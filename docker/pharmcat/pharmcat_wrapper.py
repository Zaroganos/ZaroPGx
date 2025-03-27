#!/usr/bin/env python3
"""
Simplified PharmCAT Wrapper Service for debugging
"""

import os
import json
import logging
import subprocess
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

print(f"Starting PharmCAT wrapper service with DATA_DIR={DATA_DIR}")
print(f"PharmCAT JAR location: {PHARMCAT_JAR}")

# Simple home endpoint
@app.route('/', methods=['GET'])
def home():
    return jsonify({
        "service": "PharmCAT Wrapper",
        "status": "running",
        "endpoints": ["/", "/health", "/genotype"]
    })

# Health check endpoint
@app.route('/health', methods=['GET'])
def health_check():
    """API endpoint to check if the service is running."""
    logger.info("Health check called")
    
    # Check if the JAR file exists
    jar_exists = os.path.exists(PHARMCAT_JAR)
    
    return jsonify({
        "status": "ok",
        "service": "pharmcat-wrapper",
        "java_version": subprocess.check_output(["java", "-version"], stderr=subprocess.STDOUT).decode(),
        "pharmcat_jar_exists": jar_exists,
        "pharmcat_jar_path": PHARMCAT_JAR,
        "data_dir_exists": os.path.exists(DATA_DIR),
        "data_dir_contents": os.listdir(DATA_DIR)
    })

@app.route('/genotype', methods=['POST'])
def process_genotype():
    """
    API endpoint to process VCF files with PharmCAT.
    """
    try:
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
        
        # For debugging, just return success without actually calling PharmCAT
        return jsonify({
            "success": True,
            "message": "File uploaded successfully (PharmCAT processing disabled for debugging)",
            "file_path": file_path,
            "file_size": os.path.getsize(file_path)
        })
    
    except Exception as e:
        logger.error(f"Error processing request: {str(e)}")
        return jsonify({
            "success": False,
            "message": f"Error: {str(e)}"
        })

if __name__ == '__main__':
    logger.info("Starting PharmCAT wrapper service on 0.0.0.0:5000")
    # Start Flask in debug mode
    app.run(host='0.0.0.0', port=5000, debug=True) 