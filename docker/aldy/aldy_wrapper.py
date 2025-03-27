#!/usr/bin/env python3
"""
Aldy wrapper script for CYP2D6 genotyping from VCF files.
This script provides a RESTful API endpoint to process VCF files
and return CYP2D6 genotype calls.
"""

import os
import json
import logging
import subprocess
import tempfile
from typing import Dict, Any, Optional
from flask import Flask, request, jsonify

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Initialize Flask app
app = Flask(__name__)

# Data directory for VCF files
DATA_DIR = os.environ.get("DATA_DIR", "/data")
os.makedirs(DATA_DIR, exist_ok=True)

def run_aldy(vcf_path: str, gene: str = "CYP2D6") -> Dict[str, Any]:
    """
    Run Aldy genotyping on a VCF file.
    
    Args:
        vcf_path: Path to the VCF file
        gene: Gene to genotype (default: CYP2D6)
        
    Returns:
        Dictionary containing Aldy results
    """
    try:
        logger.info(f"Running Aldy on {vcf_path} for gene {gene}")
        
        # Create a temporary file for output
        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as temp:
            output_file = temp.name
        
        # Run Aldy command
        cmd = [
            "aldy", "genotype",
            "--gene", gene,
            "--output", output_file,
            "--json",
            vcf_path
        ]
        
        process = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        logger.info(f"Aldy stdout: {process.stdout}")
        
        # Read the results JSON file
        with open(output_file, 'r') as f:
            results = json.load(f)
        
        # Clean up temporary file
        os.unlink(output_file)
        
        return {
            "gene": gene,
            "results": results,
            "diplotype": extract_diplotype(results, gene),
            "status": "success"
        }
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running Aldy: {e.stderr}")
        return {
            "gene": gene,
            "status": "error",
            "error": str(e),
            "stderr": e.stderr
        }
    except Exception as e:
        logger.error(f"Error processing Aldy results: {str(e)}")
        return {
            "gene": gene,
            "status": "error",
            "error": str(e)
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

@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint."""
    return jsonify({"status": "healthy"})

@app.route('/genotype', methods=['POST'])
def genotype():
    """
    Process a VCF file and return CYP2D6 genotype.
    
    Request should contain:
    - file: VCF file
    - gene: (optional) Gene to genotype (default: CYP2D6)
    """
    if 'file' not in request.files:
        return jsonify({"error": "No file provided"}), 400
    
    file = request.files['file']
    if not file.filename:
        return jsonify({"error": "Empty file provided"}), 400
    
    if not file.filename.endswith('.vcf'):
        return jsonify({"error": "Invalid file format. Must be VCF."}), 400
    
    gene = request.form.get('gene', 'CYP2D6')
    
    # Save the file
    file_path = os.path.join(DATA_DIR, f"{os.path.basename(file.filename)}")
    file.save(file_path)
    
    # Process with Aldy
    results = run_aldy(file_path, gene)
    
    return jsonify(results)

if __name__ == "__main__":
    # Run the Flask application
    app.run(host='0.0.0.0', port=5000) 