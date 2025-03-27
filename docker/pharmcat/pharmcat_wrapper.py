#!/usr/bin/env python3
"""
PharmCAT Wrapper Service

This script provides a Flask-based API for interacting with PharmCAT.
It handles VCF file uploads and genotype calling through PharmCAT's Java API.
"""

import os
import json
import logging
import subprocess
import tempfile
import shutil
import time
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
os.makedirs(DATA_DIR, exist_ok=True)

print("Starting PharmCAT wrapper service...")

# Make sure our health endpoint is super simple and reliable
@app.route('/health', methods=['GET'])
def health_check():
    """API endpoint to check if the service is running."""
    logger.info("Health check called")
    return jsonify({
        "status": "ok",
        "service": "pharmcat-wrapper",
        "timestamp": time.time()
    })

def run_pharmcat(vcf_path: str) -> Dict[str, Any]:
    """
    Run PharmCAT on the provided VCF file and return results.
    
    Args:
        vcf_path: Path to the VCF file
        
    Returns:
        Dictionary containing PharmCAT results or error details
    """
    try:
        # Create a temporary directory for PharmCAT output
        output_dir = tempfile.mkdtemp(dir=DATA_DIR)
        
        # Construct PharmCAT command
        cmd = [
            "java", "-jar", "/pharmcat/pharmcat.jar", "report",
            "-vcf", vcf_path,
            "-o", output_dir
        ]
        
        logger.info(f"Running PharmCAT with command: {' '.join(cmd)}")
        result = subprocess.run(
            cmd, 
            capture_output=True, 
            text=True,
            check=True
        )
        
        # Parse PharmCAT output JSON
        output_json_path = Path(output_dir) / "pharmcat_report.json"
        if output_json_path.exists():
            with open(output_json_path, 'r') as f:
                pharmcat_results = json.load(f)
                return {
                    "success": True,
                    "message": "PharmCAT analysis completed successfully",
                    "results": pharmcat_results
                }
        else:
            return {
                "success": False,
                "message": "PharmCAT completed but no results file found"
            }
    except subprocess.CalledProcessError as e:
        logger.error(f"PharmCAT execution failed: {e.stderr}")
        return {
            "success": False,
            "message": "PharmCAT execution failed",
            "error": e.stderr
        }
    except Exception as e:
        logger.error(f"Error running PharmCAT: {str(e)}")
        return {
            "success": False,
            "message": f"Error running PharmCAT: {str(e)}"
        }
    finally:
        # Clean up temporary directory
        if 'output_dir' in locals():
            shutil.rmtree(output_dir, ignore_errors=True)

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
    gene_calls = pharmcat_data.get("genotypes", [])
    
    for gene in gene_calls:
        symbol = gene.get("gene", "")
        diplotype = gene.get("diplotype", "")
        if symbol and diplotype:
            genotypes[symbol] = diplotype
            
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
    if pharmcat_results.get("success", False):
        genotypes = extract_genotypes(pharmcat_results)
        pharmcat_results["genotypes"] = genotypes
    
    return jsonify(pharmcat_results)

if __name__ == '__main__':
    logger.info("Starting PharmCAT wrapper service on 0.0.0.0:5000")
    # Add host='0.0.0.0' to make sure it's accessible from outside the container
    app.run(host='0.0.0.0', port=5000, debug=False) 