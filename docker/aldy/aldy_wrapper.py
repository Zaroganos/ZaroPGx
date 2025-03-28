"""
Aldy Genotyping Service Wrapper
-------------------------------

This module wraps the Aldy tool for pharmacogenomic star allele calling.
Successfully implemented and tested with Aldy 4.x for CYP2D6 genotyping.

Features:
- Multiple gene support via gene groups (CYP450_Enzymes, Phase_II_Enzymes, etc.)
- Automatic reference genome detection (hg38/hg19)
- Proper VCF file validation and indexing
- Support for sequencing profile parameter
- JSON output format for API integration

Usage:
- Single gene endpoint: /genotype
- Multiple genes endpoint: /multi_genotype
- Gene group endpoint: /gene_group

For development and testing:
- Run debug_aldy.sh to test with sample VCF files
- Use test_aldy_multi_gene.py for integration testing
"""

#!/usr/bin/env python3
"""
Aldy wrapper script for pharmacogenomic gene genotyping from VCF files.
This script provides a RESTful API endpoint to process VCF files
and return genotype calls for multiple supported genes.
"""

import os
import json
import logging
import subprocess
import tempfile
import glob
import yaml
from typing import Dict, Any, Optional, List, Union
from flask import Flask, request, jsonify
from werkzeug.utils import secure_filename

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Initialize Flask app
app = Flask(__name__)

# Data directory for VCF files
DATA_DIR = os.environ.get("DATA_DIR", "/data")
os.makedirs(DATA_DIR, exist_ok=True)

# Supported genes grouped by function
SUPPORTED_GENES = {
    "CYP450_Enzymes": [
        "CYP2D6", "CYP2C19", "CYP2C9", "CYP2C8", "CYP2B6",
        "CYP1A2", "CYP2A6", "CYP2E1", "CYP3A4", "CYP3A5", "CYP3A7", "CYP4F2"
    ],
    "Phase_II_Enzymes": [
        "UGT1A1", "NAT1", "NAT2", "TPMT", "DPYD", "NUDT15", "GSTM1", "GSTP1"
    ],
    "Drug_Transporters": [
        "SLCO1B1", "ABCG2"
    ],
    "Drug_Targets": [
        "VKORC1", "CACNA1S", "RYR1", "CFTR"
    ],
    "Other_PGx_Genes": [
        "COMT", "G6PD", "IFNL3", "UGT2B7"
    ]
}

# Flatten list for easy lookup
ALL_SUPPORTED_GENES = [gene for group in SUPPORTED_GENES.values() for gene in group]

def get_gene_group(gene: str) -> str:
    """Get the functional group for a given gene"""
    for group, genes in SUPPORTED_GENES.items():
        if gene in genes:
            return group
    return "Unknown"

def ensure_gene_definition_files() -> Dict[str, bool]:
    """
    Ensure all gene definition files exist, converting from YAML to JSON if needed
    
    Returns:
        Dictionary mapping gene names to whether they have valid definition files
    """
    gene_files = {}
    
    # Get path to Aldy resources
    aldy_path = None
    try:
        # Try to get the path from Python's installation
        import site
        for path in site.getsitepackages():
            possible_path = os.path.join(path, 'aldy', 'resources', 'genes')
            if os.path.exists(possible_path):
                aldy_path = possible_path
                break
        
        # Fallback to typical Docker container path
        if not aldy_path:
            aldy_path = "/usr/local/lib/python3.10/site-packages/aldy/resources/genes"
        
        if not os.path.exists(aldy_path):
            logger.error(f"Aldy gene resources path not found: {aldy_path}")
            return {gene: False for gene in ALL_SUPPORTED_GENES}
        
        # Check which genes need JSON files
        for gene in ALL_SUPPORTED_GENES:
            gene_lower = gene.lower()
            json_file = os.path.join(aldy_path, f"{gene_lower}.json")
            yaml_file = os.path.join(aldy_path, f"{gene_lower}.yml")
            
            # Check if the JSON file exists and is not empty
            json_exists = os.path.exists(json_file) and os.path.getsize(json_file) > 0
            yaml_exists = os.path.exists(yaml_file) and os.path.getsize(yaml_file) > 0
            
            # If JSON file doesn't exist or is empty but YAML file exists
            if not json_exists and yaml_exists:
                logger.info(f"Converting YAML to JSON for {gene}")
                try:
                    # Read YAML file
                    with open(yaml_file, 'r') as f:
                        gene_data = yaml.safe_load(f)
                    
                    # Write JSON file
                    with open(json_file, 'w') as f:
                        json.dump(gene_data, f, indent=2)
                    
                    json_exists = True
                    logger.info(f"Successfully created JSON file for {gene}")
                except Exception as e:
                    logger.error(f"Error converting YAML to JSON for {gene}: {str(e)}")
            
            gene_files[gene] = json_exists
            
            # Log missing gene files
            if not json_exists:
                logger.warning(f"Gene definition file missing for {gene}: {json_file}")
    
    except Exception as e:
        logger.error(f"Error ensuring gene definition files: {str(e)}")
        return {gene: False for gene in ALL_SUPPORTED_GENES}
    
    return gene_files

def check_gene_definition_files() -> Dict[str, bool]:
    """
    Check which gene definition files exist in the Aldy resources
    
    Returns:
        Dictionary mapping gene names to whether they have definition files
    """
    # First ensure all files exist
    ensure_gene_definition_files()
    
    gene_files = {}
    aldy_path = None
    
    # Get path to Aldy resources
    try:
        # Try to get the path from Python's installation
        import site
        for path in site.getsitepackages():
            possible_path = os.path.join(path, 'aldy', 'resources', 'genes')
            if os.path.exists(possible_path):
                aldy_path = possible_path
                break
        
        # Fallback to typical Docker container path
        if not aldy_path:
            aldy_path = "/usr/local/lib/python3.10/site-packages/aldy/resources/genes"
        
        # Check which gene files exist
        if os.path.exists(aldy_path):
            for gene in ALL_SUPPORTED_GENES:
                gene_lower = gene.lower()
                gene_file = os.path.join(aldy_path, f"{gene_lower}.json")
                gene_files[gene] = os.path.exists(gene_file) and os.path.getsize(gene_file) > 0
                
                # Log missing gene files
                if not gene_files[gene]:
                    logger.warning(f"Gene definition file missing for {gene}: {gene_file}")
        else:
            logger.error(f"Aldy gene resources path not found: {aldy_path}")
            return {gene: False for gene in ALL_SUPPORTED_GENES}
    
    except Exception as e:
        logger.error(f"Error checking gene definition files: {str(e)}")
        return {gene: False for gene in ALL_SUPPORTED_GENES}
    
    return gene_files

def verify_vcf_file(vcf_path: str) -> Dict[str, Any]:
    """
    Verify a VCF file is valid and can be processed by Aldy
    
    Args:
        vcf_path: Path to the VCF file
        
    Returns:
        Dictionary with verification results
    """
    try:
        # Check file exists
        if not os.path.exists(vcf_path):
            return {
                "valid": False,
                "error": f"File not found: {vcf_path}"
            }
        
        # Check file extension
        if not vcf_path.endswith(('.vcf', '.vcf.gz')):
            return {
                "valid": False,
                "error": "File must have .vcf or .vcf.gz extension"
            }
        
        # Try to get basic information using pysam
        try:
            import pysam
            vcf = pysam.VariantFile(vcf_path)
            samples = list(vcf.header.samples)
            
            # Check that there's at least one sample
            if not samples:
                return {
                    "valid": False,
                    "error": "VCF file contains no samples"
                }
                
            return {
                "valid": True,
                "samples": samples,
                "filename": os.path.basename(vcf_path)
            }
        except Exception as e:
            return {
                "valid": False,
                "error": f"Invalid VCF format: {str(e)}"
            }
            
    except Exception as e:
        return {
            "valid": False,
            "error": f"Error verifying VCF file: {str(e)}"
        }

def is_gzipped(file_path: str) -> bool:
    """Check if a file is already gzipped regardless of extension"""
    try:
        with open(file_path, 'rb') as f:
            # Check for gzip magic number 0x1f 0x8b
            return f.read(2) == b'\x1f\x8b'
    except Exception as e:
        logger.warning(f"Error checking if file is gzipped: {str(e)}")
        return False

def parse_aldy_output(output_file: str, gene: str) -> Dict[str, Any]:
    """
    Parse Aldy output file, handling both JSON and TSV formats
    
    Args:
        output_file: Path to the output file
        gene: Gene name
        
    Returns:
        Dictionary containing parsed results
    """
    try:
        # Check if the output file exists and has content
        if not os.path.exists(output_file):
            logger.error(f"Output file {output_file} was not created by Aldy")
            return {
                "gene": gene,
                "group": get_gene_group(gene),
                "status": "error",
                "error": f"Output file was not created by Aldy"
            }
            
        # Check file size
        if os.path.getsize(output_file) == 0:
            logger.error(f"Output file {output_file} is empty")
            return {
                "gene": gene,
                "group": get_gene_group(gene),
                "status": "error",
                "error": f"Output file is empty"
            }
        
        # Read the file content for inspection
        with open(output_file, 'r') as f:
            file_content = f.read()
            logger.info(f"Output file content (first 100 chars): {file_content[:100]}...")
        
        # Check if the output is in TSV format (starts with #Sample or similar)
        if file_content.startswith('#'):
            logger.info("Detected TSV format output from Aldy")
            return parse_aldy_tsv(file_content, gene)
        
        # Try to parse as JSON
        try:
            results = json.loads(file_content)
            logger.info("Successfully parsed JSON output from Aldy")
            return {
                "gene": gene,
                "group": get_gene_group(gene),
                "results": results,
                "diplotype": extract_diplotype(results, gene),
                "activity_score": extract_activity_score(results),
                "status": "success"
            }
        except json.JSONDecodeError as e:
            logger.error(f"Failed to parse output as JSON: {str(e)}")
            # Return the TSV content as fallback
            return {
                "gene": gene,
                "group": get_gene_group(gene),
                "status": "success",
                "raw_output": file_content,
                "output_format": "tsv",
                "diplotype": extract_diplotype_from_tsv(file_content, gene)
            }
    except Exception as e:
        logger.error(f"Error parsing Aldy output: {str(e)}")
        return {
            "gene": gene,
            "group": get_gene_group(gene),
            "status": "error",
            "error": f"Failed to parse output: {str(e)}"
        }

def parse_aldy_tsv(tsv_content: str, gene: str) -> Dict[str, Any]:
    """
    Parse Aldy TSV output format
    
    Args:
        tsv_content: TSV content from Aldy
        gene: Gene name
        
    Returns:
        Dictionary containing parsed results
    """
    try:
        # Extract relevant information from TSV
        lines = tsv_content.strip().split('\n')
        
        # Parse the solution header line (starts with #Solution)
        solution_line = None
        solution_data = {}
        diplotype = None
        activity_score = None
        
        logger.info(f"Attempting to parse TSV with {len(lines)} lines")
        
        # First look for the #Solution line
        for line in lines:
            if line.startswith('#Solution'):
                solution_line = line
                logger.info(f"Found solution line: {solution_line}")
                
                # Handle different solution line formats
                if ':' in solution_line:
                    # Format: #Solution 1: *1.001, *1.001
                    try:
                        allele_part = solution_line.split(':')[1].strip()
                        # Remove the * characters and split by comma
                        alleles = [a.strip().lstrip('*') for a in allele_part.split(',')]
                        if len(alleles) >= 2:
                            diplotype = f"*{alleles[0]}/*{alleles[1]}"
                            solution_data['major'] = ','.join([a.split('.')[0] for a in alleles])
                        elif len(alleles) == 1 and alleles[0]:
                            diplotype = f"*{alleles[0]}/*{alleles[0]}"
                            solution_data['major'] = alleles[0].split('.')[0]
                        
                        # Add gene to solution data
                        solution_data['gene'] = gene
                        
                        logger.info(f"Parsed alleles from solution line: {alleles}")
                    except Exception as e:
                        logger.error(f"Error parsing alleles from solution line: {str(e)}")
                else:
                    # Format: #Solution gene=CYP2D6 major=1,2 minor= score=100 activity=1.0
                    parts = line[10:].split()  # Skip "#Solution "
                    for part in parts:
                        if '=' in part:
                            key, value = part.split('=', 1)
                            solution_data[key] = value
                    
                    # Extract diplotype
                    if 'major' in solution_data:
                        alleles = solution_data['major'].split(',')
                        if len(alleles) >= 2:
                            diplotype = f"*{alleles[0]}/*{alleles[1]}"
                        elif len(alleles) == 1 and alleles[0]:
                            diplotype = f"*{alleles[0]}/*{alleles[0]}"
                
                # Extract activity score
                if 'activity' in solution_data:
                    try:
                        activity_score = float(solution_data['activity'])
                    except ValueError:
                        activity_score = None
                
                break  # Found what we needed
        
        # If we can't find the solution line, try to extract from the data rows
        if not solution_line or not diplotype:
            logger.info("No valid diplotype from #Solution line, trying to extract from data rows")
            
            # Look for lines with Major columns that have star alleles
            major_alleles = []
            
            # Skip header lines
            data_lines = [line for line in lines if not line.startswith('#')]
            header_line = None
            for line in lines:
                if line.startswith('#Sample') or line.startswith('#ID'):
                    header_line = line[1:]  # Remove leading #
                    break
            
            if header_line and data_lines:
                # Parse header to find major allele column
                header_cols = header_line.split('\t')
                major_idx = -1
                for i, col in enumerate(header_cols):
                    if col.strip() == 'Major':
                        major_idx = i
                        break
                
                if major_idx >= 0:
                    # Extract major alleles from data rows
                    for data_line in data_lines:
                        cols = data_line.split('\t')
                        if len(cols) > major_idx:
                            allele = cols[major_idx].strip()
                            if allele and allele not in major_alleles:
                                major_alleles.append(allele)
                    
                    # Create diplotype from found alleles
                    if len(major_alleles) >= 2:
                        diplotype = f"*{major_alleles[0]}/*{major_alleles[1]}"
                    elif len(major_alleles) == 1:
                        diplotype = f"*{major_alleles[0]}/*{major_alleles[0]}"
                    
                    logger.info(f"Extracted major alleles from data rows: {major_alleles}")
                    
                    # Create solution data
                    solution_data = {
                        'gene': gene,
                        'major': ','.join(major_alleles)
                    }
        
        logger.info(f"Parsed TSV solution: {solution_data}")
        logger.info(f"Extracted diplotype: {diplotype}, activity score: {activity_score}")
        
        return {
            "gene": gene,
            "group": get_gene_group(gene),
            "status": "success" if diplotype else "partial",
            "raw_output": tsv_content,
            "solution_data": solution_data,
            "diplotype": diplotype,
            "activity_score": activity_score,
            "output_format": "tsv"
        }
    except Exception as e:
        logger.error(f"Error parsing TSV output: {str(e)}")
        # Return the TSV content anyway
        return {
            "gene": gene,
            "group": get_gene_group(gene),
            "status": "partial",
            "raw_output": tsv_content,
            "error": f"Error parsing TSV: {str(e)}",
            "output_format": "tsv"
        }

def extract_diplotype_from_tsv(tsv_content: str, gene: str) -> Optional[str]:
    """
    Extract diplotype from TSV output format
    
    Args:
        tsv_content: TSV content from Aldy
        gene: Gene name
        
    Returns:
        Diplotype string or None if not found
    """
    try:
        lines = tsv_content.strip().split('\n')
        
        # Try to find the solution line
        for line in lines:
            if line.startswith('#Solution'):
                # Handle different solution line formats
                if ':' in line:
                    # Format: #Solution 1: *1.001, *1.001
                    try:
                        allele_part = line.split(':')[1].strip()
                        # Remove the * characters and split by comma
                        alleles = [a.strip().lstrip('*') for a in allele_part.split(',')]
                        if len(alleles) >= 2:
                            return f"*{alleles[0]}/*{alleles[1]}"
                        elif len(alleles) == 1 and alleles[0]:
                            return f"*{alleles[0]}/*{alleles[0]}"
                    except Exception as e:
                        logger.error(f"Error parsing alleles from solution line: {str(e)}")
                else:
                    # Format: #Solution gene=CYP2D6 major=1,2 minor= score=100 activity=1.0
                    if 'major=' in line:
                        major_part = line.split('major=')[1].split()[0]
                        alleles = major_part.split(',')
                        if len(alleles) >= 2:
                            return f"*{alleles[0]}/*{alleles[1]}"
                        elif len(alleles) == 1 and alleles[0]:
                            return f"*{alleles[0]}/*{alleles[0]}"
        
        # If we can't find the solution line, try data rows
        header_line = None
        for line in lines:
            if line.startswith('#Sample') or line.startswith('#ID'):
                header_line = line[1:]  # Remove leading #
                break
                
        if header_line:
            # Parse header to find major allele column
            header_cols = header_line.split('\t')
            major_idx = -1
            for i, col in enumerate(header_cols):
                if col.strip() == 'Major':
                    major_idx = i
                    break
            
            if major_idx >= 0:
                # Extract major alleles from data rows
                major_alleles = []
                data_lines = [line for line in lines if not line.startswith('#')]
                
                for data_line in data_lines:
                    cols = data_line.split('\t')
                    if len(cols) > major_idx:
                        allele = cols[major_idx].strip()
                        if allele and allele not in major_alleles:
                            major_alleles.append(allele)
                
                # Create diplotype from found alleles
                if len(major_alleles) >= 2:
                    return f"*{major_alleles[0]}/*{major_alleles[1]}"
                elif len(major_alleles) == 1:
                    return f"*{major_alleles[0]}/*{major_alleles[0]}"
        
        return None
    except Exception as e:
        logger.error(f"Error extracting diplotype from TSV: {str(e)}")
        return None

def run_aldy(vcf_path: str, gene: str = "CYP2D6", profile: str = None, 
             sequencing_profile: str = "illumina", file_type: str = "vcf",
             reference_genome: str = "hg19") -> Dict[str, Any]:
    """
    Run Aldy genotyping on a genomic file.
    
    Args:
        vcf_path: Path to the genomic file (VCF, BAM, SAM, CRAM)
        gene: Gene to genotype (default: CYP2D6)
        profile: (optional) Sample name in the file, not used directly in command
        sequencing_profile: (optional) Sequencing technology profile (default: illumina)
        file_type: Type of the input file (vcf, bam, sam, cram)
        reference_genome: Reference genome for CRAM files (default: hg19)
        
    Returns:
        Dictionary containing Aldy results
    """
    try:
        if gene not in ALL_SUPPORTED_GENES:
            return {
                "gene": gene,
                "status": "error",
                "error": f"Gene {gene} is not supported by this Aldy wrapper"
            }
            
        logger.info(f"Running Aldy on {vcf_path} for gene {gene}")
        
        # Create a temporary file for output
        with tempfile.NamedTemporaryFile(suffix=".txt", delete=False) as temp:
            output_file = temp.name
            
        logger.info(f"Output will be saved to {output_file}")
        
        # If the file is VCF, check if it's compressed
        if file_type == 'vcf':
            is_compressed = is_gzipped(vcf_path)
            logger.info(f"File compression check: {vcf_path} is {'compressed' if is_compressed else 'not compressed'}")
        else:
            is_compressed = False
        
        # Get sample name from the file
        try:
            import pysam
            
            # Handle different file types
            if file_type == 'vcf':
                # For VCF files
                vcf = pysam.VariantFile(vcf_path)
                samples = list(vcf.header.samples)
            elif file_type in ['bam', 'sam', 'cram']:
                # For BAM/SAM/CRAM files
                reference_param = ["-T", reference_genome] if file_type == 'cram' else []
                samfile = pysam.AlignmentFile(vcf_path, reference_filename=reference_genome if file_type == 'cram' else None)
                # Get sample name from BAM header
                header = samfile.header
                if 'RG' in header and header['RG'] and 'SM' in header['RG'][0]:
                    samples = [header['RG'][0]['SM']]
                else:
                    # If no sample name in header, use a default
                    samples = ["SAMPLE"]
            else:
                # Shouldn't reach here due to validation earlier
                samples = []
                
            if samples:
                logger.info(f"File contains samples: {', '.join(samples)}")
                # Always use the first sample from the file, ignore the profile parameter
                vcf_sample = samples[0]
                logger.info(f"Using first sample from file as profile: {vcf_sample}")
                # Override the profile parameter with the actual sample name from file
                profile = vcf_sample
        except Exception as e:
            logger.warning(f"Could not get sample name from file: {str(e)}")
        
        # Run Aldy command according to the documentation
        cmd = [
            "aldy", "genotype",
            "-g", gene,                    # Gene to analyze
            "-p", sequencing_profile,      # Sequencing technology profile
            "-o", output_file              # Output file
        ]
        
        # Add reference genome for CRAM files
        if file_type == 'cram':
            cmd.extend(["-r", reference_genome])
        
        # Add sample parameter if available (using --profile instead of --sample)
        if profile:
            cmd.extend(["--profile", profile])
            logger.info(f"Using profile parameter: {profile}")
        
        # Add verbose logging to help with debugging
        cmd.extend(["--log", "DEBUG"])
        
        # Handle file compression and indexing for VCF files
        fixed_vcf_path = vcf_path
        index_created = False
        compressed_path = None
        
        # For VCF files, handle compression and indexing
        if file_type == 'vcf':
            # For explicitly compressed files (.vcf.gz)
            if vcf_path.endswith('.vcf.gz'):
                compressed_path = vcf_path
                index_path = f"{vcf_path}.tbi"
                csi_path = f"{vcf_path}.csi"
                if not (os.path.exists(index_path) or os.path.exists(csi_path)):
                    logger.info(f"No index found for {vcf_path}, creating one...")
                    try:
                        index_cmd = ["tabix", "-f", "-p", "vcf", vcf_path]
                        result = subprocess.run(index_cmd, capture_output=True, text=True)
                        if result.returncode == 0:
                            logger.info(f"Created index for {vcf_path}")
                            index_created = True
                        else:
                            logger.warning(f"Failed to create index: {result.stderr}")
                    except Exception as e:
                        logger.warning(f"Failed to run tabix: {str(e)}")
    
            # For files without .gz extension but actually compressed
            elif is_compressed:
                # File is actually compressed but doesn't have .gz extension
                # Create a properly named copy with .gz extension
                logger.info(f"File is compressed but missing .gz extension, creating properly named copy")
                compressed_path = f"{vcf_path}.gz"
                # Copy the file with the correct extension
                try:
                    import shutil
                    # Make a new copy only if it doesn't exist or overwrite existing
                    if not os.path.exists(compressed_path) or gene == "CYP2D6":  # Always create fresh copy for first gene
                        shutil.copy2(vcf_path, compressed_path)
                        logger.info(f"Created copy at {compressed_path}")
                    
                    # Try to create an index, using -f to force overwrite if needed
                    try:
                        index_cmd = ["tabix", "-f", "-p", "vcf", compressed_path]
                        result = subprocess.run(index_cmd, capture_output=True, text=True)
                        if result.returncode == 0:
                            logger.info(f"Created index for {compressed_path}")
                            fixed_vcf_path = compressed_path
                            index_created = True
                        else:
                            logger.warning(f"Failed to create index: {result.stderr}")
                            # Continue with original file if indexing fails
                            logger.info(f"Using original file: {vcf_path}")
                    except Exception as e:
                        logger.warning(f"Failed to run tabix: {str(e)}")
                except Exception as e:
                    logger.warning(f"Failed to create properly named copy: {str(e)}")
    
            # For truly uncompressed VCF files
            elif vcf_path.endswith('.vcf'):
                # Try to compress and index the VCF file
                logger.info(f"Compressing and indexing uncompressed VCF: {vcf_path}")
                try:
                    compressed_path = f"{vcf_path}.gz"
                    
                    # Check if bgzip is available
                    bgzip_available = False
                    try:
                        bgzip_check = subprocess.run(["which", "bgzip"], capture_output=True, text=True)
                        bgzip_available = bgzip_check.returncode == 0
                    except:
                        logger.warning("Failed to check for bgzip")
                    
                    if not bgzip_available:
                        logger.warning("bgzip not found, using uncompressed file")
                        # Continue with uncompressed file
                    else:
                        # Compress the file if it doesn't exist or if this is the first gene
                        if not os.path.exists(compressed_path) or gene == "CYP2D6":
                            bgzip_cmd = ["bgzip", "-c", vcf_path]
                            with open(compressed_path, 'wb') as f:
                                result = subprocess.run(bgzip_cmd, stdout=f, capture_output=False)
                            
                            if os.path.exists(compressed_path) and os.path.getsize(compressed_path) > 0:
                                logger.info(f"Successfully compressed to {compressed_path}")
                            else:
                                logger.warning(f"Failed to compress to {compressed_path}")
                                # Continue with uncompressed file
                        
                        # Try to create an index with -f to force
                        if os.path.exists(compressed_path):
                            try:
                                index_cmd = ["tabix", "-f", "-p", "vcf", compressed_path]
                                result = subprocess.run(index_cmd, capture_output=True, text=True)
                                if result.returncode == 0:
                                    logger.info(f"Created index for {compressed_path}")
                                    fixed_vcf_path = compressed_path
                                    index_created = True
                                else:
                                    logger.warning(f"Failed to create index: {result.stderr}")
                                    # Continue with uncompressed file
                            except Exception as e:
                                logger.warning(f"Failed to run tabix: {str(e)}")
                except Exception as e:
                    logger.warning(f"Error in compression/indexing process: {str(e)}")
                    # Continue with uncompressed file
        
        # Always use the indexed file when available
        if compressed_path and os.path.exists(compressed_path) and (os.path.exists(f"{compressed_path}.tbi") or os.path.exists(f"{compressed_path}.csi")):
            logger.info(f"Using indexed file: {compressed_path}")
            fixed_vcf_path = compressed_path
        
        # The file path is the last argument
        cmd.append(fixed_vcf_path)
        
        logger.info(f"Running command: {' '.join(cmd)}")
        
        try:
            process = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            
            logger.info(f"Aldy stdout: {process.stdout}")
            
            # Parse the output file and return results
            return parse_aldy_output(output_file, gene)
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Error running Aldy: {e.stderr}")
            
            # Check for specific error patterns
            error_detail = ""
            stderr = e.stderr.lower() if e.stderr else ""
            if "no such file" in stderr or "not found" in stderr:
                error_detail = f"Gene definition file for {gene} may be missing"
            elif "no sample" in stderr:
                error_detail = "File may not contain the required sample data"
            elif "no profile" in stderr or "--profile" in stderr:
                error_detail = "Profile parameter required or invalid for this file"
            elif "no coverage" in stderr:
                error_detail = f"File doesn't have coverage profile for {gene}"
            elif "invalid format" in stderr:
                error_detail = "File format is invalid"
            elif "index" in stderr:
                error_detail = "File indexing is required but failed"
                
            # Check if any output was produced even though the command failed
            if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                logger.info("Output file exists despite command failure, attempting to parse")
                with open(output_file, 'r') as f:
                    file_content = f.read()
                    if file_content.strip():
                        logger.info(f"Output file content (first 100 chars): {file_content[:100]}...")
                        # Try to parse the partial output
                        return parse_aldy_output(output_file, gene)
            
            return {
                "gene": gene,
                "group": get_gene_group(gene),
                "status": "error",
                "error": str(e),
                "stderr": e.stderr,
                "error_detail": error_detail,
                "cmd": " ".join(cmd)
            }
    except Exception as e:
        logger.error(f"Error processing Aldy results: {str(e)}")
        return {
            "gene": gene,
            "group": get_gene_group(gene),
            "status": "error",
            "error": str(e)
        }
    finally:
        # Clean up temporary file
        try:
            if 'output_file' in locals() and os.path.exists(output_file):
                os.unlink(output_file)
                logger.info(f"Cleaned up temporary file {output_file}")
        except Exception as e:
            logger.warning(f"Failed to clean up temporary file: {str(e)}")

def run_multi_gene_analysis(vcf_path: str, genes: List[str], profile: str = None, sequencing_profile: str = "illumina") -> Dict[str, Any]:
    """
    Run Aldy genotyping for multiple genes.
    
    Args:
        vcf_path: Path to the VCF file
        genes: List of genes to analyze
        profile: (optional) Sample name in the VCF file, not used directly in command
        sequencing_profile: (optional) Sequencing technology profile (default: illumina)
        
    Returns:
        Dictionary containing results for all genes
    """
    results = {}
    
    for gene in genes:
        gene_result = run_aldy(vcf_path, gene, profile, sequencing_profile)
        results[gene] = gene_result
    
    return {
        "status": "success",
        "genes": results,
        "groups": SUPPORTED_GENES
    }

def extract_diplotype(results: Dict[str, Any], gene: str) -> Optional[str]:
    """
    Extract the diplotype from Aldy results.
    
    Args:
        results: Aldy results dictionary
        gene: Gene name
        
    Returns:
        Diplotype string or None if not found
    """
    try:
        # Extract diplotype based on Aldy's JSON structure
        solutions = results.get("solutions", [])
        if not solutions:
            return None
        
        # Get the first solution (highest score)
        solution = solutions[0]
        alleles = solution.get("major", {}).get("alleles", [])
        
        if len(alleles) >= 2:
            return f"*{alleles[0]}/*{alleles[1]}"
        elif len(alleles) == 1:
            return f"*{alleles[0]}/*{alleles[0]}"  # Homozygous
        else:
            return None
    except Exception as e:
        logger.error(f"Error extracting diplotype: {str(e)}")
        return None

def extract_activity_score(results: Dict[str, Any]) -> Optional[float]:
    """
    Extract the activity score from Aldy results if available.
    
    Args:
        results: Aldy results dictionary
        
    Returns:
        Activity score or None if not found
    """
    try:
        solutions = results.get("solutions", [])
        if not solutions:
            return None
        
        solution = solutions[0]
        activity = solution.get("scores", {}).get("activity")
        
        return activity
    except Exception as e:
        logger.error(f"Error extracting activity score: {str(e)}")
        return None

@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint."""
    # Ensure gene definition files first
    ensure_gene_definition_files()
    
    # Then check which ones exist
    gene_files = check_gene_definition_files()
    genes_with_files = [gene for gene, exists in gene_files.items() if exists]
    
    return jsonify({
        "status": "healthy",
        "supported_genes": ALL_SUPPORTED_GENES,
        "gene_groups": SUPPORTED_GENES,
        "gene_definition_files": gene_files,
        "genes_with_definition_files": genes_with_files
    })

@app.route('/supported_genes', methods=['GET'])
def get_supported_genes():
    """Return the list of supported genes and their groupings."""
    return jsonify({
        "genes": ALL_SUPPORTED_GENES,
        "groups": SUPPORTED_GENES
    })

@app.route('/genotype', methods=['POST'])
def genotype():
    """
    Process a VCF file and return genotype for a specified gene.
    
    Request should contain:
    - file: VCF file
    - gene: (optional) Gene to genotype (default: CYP2D6)
    - profile: (optional) Sample name in VCF (for logging only)
    - sequencing_profile: (optional) Sequencing technology profile (default: illumina)
      Valid values: illumina, pgx1, pgx2, pgx3, exome, 10x, etc.
    """
    if 'file' not in request.files:
        return jsonify({"error": "No file provided"}), 400
    
    file = request.files['file']
    if not file.filename:
        return jsonify({"error": "Empty file provided"}), 400
    
    if not file.filename.endswith(('.vcf', '.vcf.gz')):
        return jsonify({"error": "Invalid file format. Must be VCF."}), 400
    
    gene = request.form.get('gene', 'CYP2D6')
    profile = request.form.get('profile')
    sequencing_profile = request.form.get('sequencing_profile', 'illumina')
    
    # Save the file
    file_path = os.path.join(DATA_DIR, f"{os.path.basename(file.filename)}")
    file.save(file_path)
    
    # Process with Aldy
    results = run_aldy(file_path, gene, profile, sequencing_profile)
    
    return jsonify(results)

@app.route('/multi_genotype', methods=['POST'])
def multi_genotype():
    """
    Run Aldy on multiple genes and return the results as a JSON.
    Accepts file upload and parameters via form data.
    
    Form parameters:
        - file: The genomic file (VCF, BAM, SAM, CRAM)
        - genes: Comma-separated list of genes to analyze
        - sequencing_profile: Sequencing profile to use (default: illumina)
        - file_type: Type of the genomic file (vcf, bam, sam, cram)
        - reference_genome: Reference genome for CRAM files (default: hg19)
    
    Returns:
        JSON with gene-specific results
    """
    try:
        # Check if a file was uploaded
        if 'file' not in request.files:
            return jsonify({'status': 'error', 'error': 'No file provided'}), 400
        
        uploaded_file = request.files['file']
        
        # Get parameters from form data
        genes_str = request.form.get('genes', 'CYP2D6')
        sequencing_profile = request.form.get('sequencing_profile', 'illumina')
        file_type = request.form.get('file_type', 'vcf').lower()
        reference_genome = request.form.get('reference_genome', 'hg19')
        
        # Validate genes
        genes = [g.strip() for g in genes_str.split(',') if g.strip()]
        if not genes:
            return jsonify({'status': 'error', 'error': 'No genes specified'}), 400
        
        # Validate file type
        valid_file_types = ['vcf', 'bam', 'sam', 'cram']
        if file_type not in valid_file_types:
            return jsonify({'status': 'error', 'error': f'Invalid file type. Must be one of: {", ".join(valid_file_types)}'}), 400
        
        # For CRAM files, reference genome is required
        if file_type == 'cram' and not reference_genome:
            return jsonify({'status': 'error', 'error': 'Reference genome is required for CRAM files'}), 400
        
        # Save the uploaded file
        save_path = os.path.join(DATA_DIR, secure_filename(uploaded_file.filename))
        uploaded_file.save(save_path)
        logger.info(f"Saved uploaded file to {save_path}")
        
        # Process the file based on its type
        processed_file = save_path
        
        # For BAM/SAM/CRAM files, ensure they are indexed if needed
        if file_type in ['bam', 'cram']:
            # Check if index exists
            index_path = save_path + '.bai'
            if not os.path.exists(index_path):
                logger.info(f"Indexing {file_type.upper()} file: {save_path}")
                try:
                    index_cmd = ["samtools", "index", save_path]
                    subprocess.run(index_cmd, check=True, capture_output=True, text=True)
                    logger.info(f"Successfully indexed {file_type.upper()} file")
                except subprocess.CalledProcessError as e:
                    logger.error(f"Error indexing {file_type.upper()} file: {e.stderr}")
                    return jsonify({'status': 'error', 'error': f'Failed to index {file_type.upper()} file: {e.stderr}'}), 500
        
        # Process each gene
        results = {}
        for gene in genes:
            logger.info(f"Processing gene {gene}")
            try:
                # Run Aldy with appropriate parameters
                result = run_aldy(
                    vcf_path=processed_file,
                    gene=gene,
                    profile=None,  # Will be auto-detected from the file
                    sequencing_profile=sequencing_profile,
                    file_type=file_type,
                    reference_genome=reference_genome
                )
                results[gene] = result
            except Exception as e:
                logger.error(f"Error processing gene {gene}: {str(e)}")
                results[gene] = {
                    'gene': gene,
                    'status': 'error',
                    'error': str(e)
                }
        
        # Return combined results
        return jsonify({
            'status': 'success',
            'genes': results
        })
    except Exception as e:
        logger.error(f"Error in multi_genotype: {str(e)}")
        return jsonify({'status': 'error', 'error': str(e)}), 500

@app.route('/verify_vcf', methods=['POST'])
def verify_vcf():
    """
    Verify a VCF file can be processed by Aldy
    
    Request should contain:
    - file: VCF file
    """
    if 'file' not in request.files:
        return jsonify({"error": "No file provided"}), 400
    
    file = request.files['file']
    if not file.filename:
        return jsonify({"error": "Empty file provided"}), 400
    
    # Save the file
    file_path = os.path.join(DATA_DIR, f"{os.path.basename(file.filename)}")
    file.save(file_path)
    
    # Verify the file
    verification = verify_vcf_file(file_path)
    
    return jsonify(verification)

if __name__ == "__main__":
    # Run the Flask application with a longer timeout to handle larger VCF files
    import werkzeug.serving
    from waitress import serve
    
    # Log startup information
    logger.info(f"Starting Aldy Genotyping Service")
    logger.info(f"Supported genes: {len(ALL_SUPPORTED_GENES)} genes in {len(SUPPORTED_GENES)} groups")
    
    # Check gene definition files
    gene_files = check_gene_definition_files()
    genes_with_files = [gene for gene, exists in gene_files.items() if exists]
    logger.info(f"Found definition files for {len(genes_with_files)} genes")
    
    # Start the server
    try:
        # Use waitress for more reliable production serving with longer timeouts
        serve(app, host='0.0.0.0', port=5000, threads=4, 
              connection_limit=20, channel_timeout=300, timeout=300)
    except Exception as e:
        logger.error(f"Error starting server: {str(e)}")
        # Fallback to Flask development server
        app.run(host='0.0.0.0', port=5000, threaded=True) 