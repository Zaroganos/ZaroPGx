#!/usr/bin/env python3
"""
PharmCAT Wrapper Service

This script provides a Flask-based API for interacting with PharmCAT.
It handles VCF file uploads and executes PharmCAT directly via the JAR file.
"""

import os
import json
import logging
import subprocess
import tempfile
import shutil
import time
import requests
from typing import Dict, Any, Optional, List
from pathlib import Path
from flask import Flask, request, jsonify

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("pharmcat_wrapper")

# Initialize Flask app
app = Flask(__name__)

# Data directory for VCF files
DATA_DIR = os.environ.get("DATA_DIR", "/data")
# Path to the PharmCAT JAR file (mounted from the PharmCAT container)
PHARMCAT_JAR = os.environ.get("PHARMCAT_JAR", "/pharmcat/pharmcat.jar")

os.makedirs(DATA_DIR, exist_ok=True)

print(f"Starting PharmCAT wrapper service - will use JAR file at {PHARMCAT_JAR}")

# Make sure our health endpoint is super simple and reliable
@app.route('/health', methods=['GET'])
def health_check():
    """API endpoint to check if the service is running."""
    logger.info("Health check called")
    # Try to verify PharmCAT availability
    try:
        # Check if the JAR file exists
        if os.path.exists(PHARMCAT_JAR):
            pharmcat_status = "JAR file available"
        else:
            pharmcat_status = f"JAR file not found at {PHARMCAT_JAR}"
    except Exception as e:
        pharmcat_status = f"Error: {str(e)}"
        
    return jsonify({
        "status": "ok",
        "service": "pharmcat-wrapper",
        "pharmcat_jar": pharmcat_status,
        "timestamp": time.time()
    })

def run_pharmcat(vcf_path: str) -> Dict[str, Any]:
    """
    Run PharmCAT on the provided VCF file by directly calling the JAR.
    
    Args:
        vcf_path: Path to the VCF file
        
    Returns:
        Dictionary containing PharmCAT results or error details
    """
    try:
        # Get just the filename for better output naming
        filename = os.path.basename(vcf_path)
        base_filename = os.path.splitext(filename)[0]
        
        # Set output directory to our data dir
        out_dir = DATA_DIR
        
        logger.info(f"Running PharmCAT on VCF: {vcf_path}")
        
        # Execute PharmCAT directly
        cmd = [
            "java", "-jar", PHARMCAT_JAR,
            "-vcf", vcf_path,
            "-o", out_dir,
            "-bf", base_filename
        ]
        
        logger.info(f"Running command: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            timeout=120
        )
        
        logger.info(f"PharmCAT execution successful: {result.stdout}")
        
        # Try loading the phenotype JSON which is likely what we want
        phenotype_json = Path(f"{DATA_DIR}/{base_filename}.phenotype.json")
        if phenotype_json.exists():
            with open(phenotype_json, 'r') as f:
                try:
                    pharmcat_results = json.load(f)
                    return {
                        "success": True,
                        "message": "PharmCAT analysis completed successfully",
                        "results": pharmcat_results
                    }
                except json.JSONDecodeError:
                    logger.error(f"Failed to parse PharmCAT phenotype JSON: {phenotype_json}")
        
        # Try loading match.json as a fallback
        match_json = Path(f"{DATA_DIR}/{base_filename}.match.json")
        if match_json.exists():
            with open(match_json, 'r') as f:
                try:
                    pharmcat_results = json.load(f)
                    return {
                        "success": True,
                        "message": "PharmCAT analysis completed successfully (no phenotype data)",
                        "results": pharmcat_results
                    }
                except json.JSONDecodeError:
                    logger.error(f"Failed to parse PharmCAT match JSON: {match_json}")
        
        # Try loading the HTML report
        html_report = Path(f"{DATA_DIR}/{base_filename}.report.html")
        if html_report.exists():
            return {
                "success": True,
                "message": "PharmCAT completed successfully (HTML report only)",
                "report_path": str(html_report),
                "stdout": result.stdout
            }
            
        # No valid output found
        return {
            "success": False,
            "message": "PharmCAT completed but no results file found",
            "stdout": result.stdout,
            "stderr": result.stderr
        }
            
    except subprocess.CalledProcessError as e:
        logger.error(f"PharmCAT execution failed: {e.stderr}")
        return {
            "success": False,
            "message": "PharmCAT execution failed",
            "error": e.stderr
        }
    except subprocess.TimeoutExpired:
        logger.error("PharmCAT execution timed out after 120 seconds")
        return {
            "success": False,
            "message": "PharmCAT execution timed out after 120 seconds"
        }
    except Exception as e:
        logger.error(f"Error running PharmCAT: {str(e)}")
        return {
            "success": False,
            "message": f"Error running PharmCAT: {str(e)}"
        }

def extract_genotypes(results: Dict[str, Any]) -> Dict[str, str]:
    """
    Extract genotype calls from PharmCAT results.
    
    Args:
        results: PharmCAT results dictionary
        
    Returns:
        Dictionary mapping gene symbols to diplotype calls
    """
    genotypes = {}
    
    if not results.get("success", False):
        return genotypes
        
    pharmcat_data = results.get("results", {})
    
    # The structure may depend on which JSON file we loaded
    # First try the standard format
    gene_calls = pharmcat_data.get("genotypes", [])
    
    if gene_calls:
        for gene in gene_calls:
            symbol = gene.get("gene", "")
            diplotype = gene.get("diplotype", "")
            if symbol and diplotype:
                genotypes[symbol] = diplotype
    else:
        # Try alternate formats that might be in the PharmCAT output
        for gene_key, gene_data in pharmcat_data.items():
            if isinstance(gene_data, dict) and "diplotype" in gene_data:
                genotypes[gene_key] = gene_data["diplotype"]
                
    return genotypes

@app.route('/genotype', methods=['POST'])
def process_genotype():
    """
    API endpoint to process VCF files with PharmCAT.
    
    Expects a POST request with a VCF file attached.
    Returns PharmCAT genotype results.
    """
    # Check if file is in request
    if 'file' not in request.files:
        return jsonify({
            "success": False,
            "message": "No file provided"
        }), 400
        
    vcf_file = request.files['file']
    
    # Validate file
    if vcf_file.filename == '':
        return jsonify({
            "success": False,
            "message": "Empty filename"
        }), 400
        
    if not vcf_file.filename.endswith(('.vcf', '.vcf.gz')):
        return jsonify({
            "success": False,
            "message": "File must be a VCF (.vcf or .vcf.gz)"
        }), 400
    
    # Save file to data directory
    file_path = os.path.join(DATA_DIR, vcf_file.filename)
    vcf_file.save(file_path)
    
    logger.info(f"VCF file saved to {file_path}")
    
    # Process file with PharmCAT
    pharmcat_results = run_pharmcat(file_path)
    
    # Add genotype summary if successful
    if pharmcat_results.get("success", False) and "results" in pharmcat_results:
        genotypes = extract_genotypes(pharmcat_results)
        pharmcat_results["genotypes"] = genotypes
    
    return jsonify(pharmcat_results)

if __name__ == '__main__':
    logger.info("Starting PharmCAT wrapper service on 0.0.0.0:5000")
    # Add host='0.0.0.0' to make sure it's accessible from outside the container
    app.run(host='0.0.0.0', port=5000, debug=False) 