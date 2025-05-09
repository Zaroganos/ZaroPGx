#!/usr/bin/env python3
"""
Simplified PharmCAT Wrapper Service for debugging
"""

import os
import json
import logging
import subprocess
import gzip
import shutil
from pathlib import Path
from flask import Flask, request, jsonify
import psutil
import time
from datetime import datetime
import re
import tempfile
import uuid
import traceback
from typing import Dict, Any, Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("pharmcat")

# Initialize Flask app
app = Flask(__name__)

# Data directory for VCF files
DATA_DIR = os.environ.get("DATA_DIR", "/data")
# Path to the PharmCAT JAR file (for backward compatibility)
PHARMCAT_JAR = os.environ.get("PHARMCAT_JAR", "/pharmcat/pharmcat.jar")
# Path to the PharmCAT pipeline directory
PHARMCAT_PIPELINE_DIR = os.environ.get("PHARMCAT_PIPELINE_DIR", "/pharmcat/pipeline")

os.makedirs(DATA_DIR, exist_ok=True)

print(f"Starting PharmCAT wrapper service with DATA_DIR={DATA_DIR}")
print(f"PharmCAT JAR location: {PHARMCAT_JAR}")

# Add these global variables after the existing ones
processing_status = {
    "current_file": None,
    "start_time": None,
    "status": "idle",
    "progress": 0,
    "last_error": None
}

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
        "service": "pharmcat",
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

@app.route('/status', methods=['GET'])
def get_status():
    """API endpoint to get the current processing status."""
    try:
        # Get container stats
        container_stats = {
            "cpu_percent": psutil.cpu_percent(),
            "memory_percent": psutil.virtual_memory().percent,
            "disk_usage": psutil.disk_usage('/').percent
        }
        
        # Get process stats if we're processing a file
        process_stats = None
        if processing_status["current_file"]:
            process_stats = {
                "file": processing_status["current_file"],
                "start_time": processing_status["start_time"].isoformat() if processing_status["start_time"] else None,
                "elapsed_time": (datetime.now() - processing_status["start_time"]).total_seconds() if processing_status["start_time"] else 0
            }
        
        return jsonify({
            "status": "ok",
            "service": "pharmcat",
            "processing_status": processing_status,
            "process_stats": process_stats,
            "container_stats": container_stats,
            "pharmcat_jar_exists": os.path.exists(PHARMCAT_JAR),
            "data_dir_contents": os.listdir(DATA_DIR)
        })
    except Exception as e:
        logger.error(f"Error getting status: {str(e)}")
        return jsonify({
            "status": "error",
            "message": str(e)
        }), 500

@app.route('/process', methods=['POST'])
def process_file():
    """
    Process uploaded VCF file with PharmCAT
    
    Returns:
        JSON response with results or error message
    """
    processing_status = {
        "status": "processing",
        "progress": 0,
        "last_error": None
    }
    
    try:
        # Get the uploaded file from the request
        if 'file' not in request.files:
            return jsonify({
                "success": False,
                "message": "No file provided"
            }), 400
            
        vcf_file = request.files['file']
        
        # Create a temporary directory for this job
        with tempfile.TemporaryDirectory() as temp_dir:
            # Get or generate a unique name for this job
            # Use reportId from request data if provided (for consistent directory naming)
            report_id = request.form.get('reportId')
            if report_id:
                base_name = report_id
                logger.info(f"Using provided report ID for base name: {base_name}")
            else:
                # Generate a random ID if none provided
                base_name = str(uuid.uuid4())[:8]
                logger.info(f"Generated random base name: {base_name}")
            
            # Save the uploaded file
            vcf_path = os.path.join(temp_dir, f"{base_name}.vcf")
            vcf_file.save(vcf_path)
            logger.info(f"Saved uploaded file to {vcf_path}")
            
            # Update processing status
            processing_status.update({
                "status": "processing",
                "progress": 10,
                "message": "File uploaded, starting PharmCAT pipeline"
            })
            
            try:
                # Run PharmCAT
                logger.info(f"Running PharmCAT on {vcf_path}")
                
                # Configure paths for PharmCAT
                output_dir = temp_dir
                
                # Before running PharmCAT, set up Java options for memory management
                java_options = "-Xmx2g"
                
                # Build the PharmCAT command using pharmcat_pipeline as per official documentation
                # Format: pharmcat_pipeline [input file] -o [output dir] [other options]
                pharmcat_cmd = [
                    "pharmcat_pipeline",
                    "-G",  # Bypass gVCF check
                    "-o", output_dir,  # Output directory
                    "-v",  # Verbose output
                    "-reporterHtml",    # Generate HTML report
                    "-reporterJson",    # Generate JSON report 
                    "-reporterCallsOnlyTsv",  # Generate TSV report with only calls data
                    vcf_path  # Input file should be the last argument
                ]
                
                # Set environment variables
                env = os.environ.copy()
                env["JAVA_TOOL_OPTIONS"] = "-Xmx4g -XX:+UseG1GC"
                env["PHARMCAT_LOG_LEVEL"] = "DEBUG"
                
                # Log important environment info for debugging
                logger.info(f"PHARMCAT_JAR location: {PHARMCAT_JAR}")
                logger.info(f"PHARMCAT_PIPELINE_DIR: {PHARMCAT_PIPELINE_DIR}")
                logger.info(f"PATH environment: {env.get('PATH', 'Not set')}")
                
                # Check if pharmcat_pipeline exists and is executable
                try:
                    subprocess.run(["which", "pharmcat_pipeline"], check=True, capture_output=True, text=True)
                    logger.info("pharmcat_pipeline command found in PATH")
                except subprocess.CalledProcessError:
                    logger.warning("pharmcat_pipeline command NOT found in PATH")
                
                # Run PharmCAT pipeline
                logger.info(f"Executing PharmCAT command: {' '.join(pharmcat_cmd)}")
                process = subprocess.run(
                    pharmcat_cmd,
                    check=True,
                    capture_output=True,
                    text=True,
                    env=env,
                    timeout=300  # 5 minute timeout
                )
                
                logger.info(f"PharmCAT process completed with output: {process.stdout}")
                
                # Check for results files with the correct base name
                # Look for all report formats
                actual_files = [f for f in os.listdir(temp_dir) if f.endswith(('.json', '.html', '.tsv'))]
                logger.info(f"PharmCAT output files found: {actual_files}")
                
                # Find all report files by format
                report_json_file = next((f for f in actual_files if f.endswith('.report.json')), None)
                report_html_file = next((f for f in actual_files if f.endswith('.report.html')), None)
                report_tsv_file = next((f for f in actual_files if f.endswith('.report.tsv')), None)
                
                if not report_json_file:
                    error_msg = "Required PharmCAT report.json file not found"
                    logger.error(error_msg)
                    return jsonify({
                        "success": False,
                        "message": error_msg,
                        "files_found": actual_files
                    }), 500
                
                # Create report directory if it doesn't exist
                reports_dir = Path("/data/reports")
                reports_dir.mkdir(parents=True, exist_ok=True)
                
                # Create a patient-specific directory for this job
                patient_dir = reports_dir / base_name
                patient_dir.mkdir(parents=True, exist_ok=True)
                logger.info(f"Created patient-specific directory: {patient_dir}")
                
                # Copy PharmCAT HTML report to patient directory if available
                if report_html_file:
                    dest_html_path = patient_dir / f"{base_name}_pgx_report.html"
                    src_html_path = Path(temp_dir) / report_html_file
                    
                    logger.info(f"Copying PharmCAT report from {src_html_path} to {dest_html_path}")
                    
                    if os.path.exists(src_html_path):
                        shutil.copy2(src_html_path, dest_html_path)
                        logger.info(f"HTML report copied to {dest_html_path}")
                
                # Copy PharmCAT TSV report to /data/reports if available
                if report_tsv_file:
                    dest_tsv_path = patient_dir / f"{base_name}_pgx_report.tsv"
                    src_tsv_path = Path(temp_dir) / report_tsv_file
                    
                    logger.info(f"Copying PharmCAT TSV report from {src_tsv_path} to {dest_tsv_path}")
                    
                    if os.path.exists(src_tsv_path):
                        shutil.copy2(src_tsv_path, dest_tsv_path)
                        logger.info(f"TSV report copied to {dest_tsv_path}")
                
                # Copy the report.json file to reports directory for inspection
                if report_json_file:
                    dest_json_path = patient_dir / f"{base_name}_pgx_report.json"
                    src_json_path = Path(temp_dir) / report_json_file
                    
                    logger.info(f"Copying JSON report from {src_json_path} to {dest_json_path}")
                    
                    # Check if the JSON file actually exists before trying to copy it
                    if os.path.exists(src_json_path):
                        try:
                            shutil.copy2(src_json_path, dest_json_path)
                            logger.info(f"JSON report copied to {dest_json_path}")
                            # Verify the file was copied successfully
                            if os.path.exists(dest_json_path):
                                logger.info(f"Verified JSON report exists at {dest_json_path}")
                            else:
                                logger.error(f"Failed to copy JSON report to {dest_json_path}")
                        except Exception as e:
                            logger.error(f"Error copying JSON report: {str(e)}")
                    else:
                        logger.error(f"JSON report not found at {src_json_path}")
                        # List all files in the temp directory for debugging
                        all_files = os.listdir(temp_dir)
                        logger.info(f"All files in temp directory: {all_files}")
                        
                        # Create a backup JSON file with basic structure for inspection
                        logger.info("Creating a backup JSON file for inspection")
                        backup_json_path = patient_dir / f"{base_name}_pgx_report.json"
                        try:
                            # Get data from the phenotype.json file if available
                            phenotype_file = next((f for f in actual_files if f.endswith('.phenotype.json')), None)
                            if phenotype_file:
                                logger.info(f"Using phenotype file {phenotype_file} for backup")
                                with open(Path(temp_dir) / phenotype_file, 'r') as f:
                                    phenotype_data = json.load(f)
                                    
                                # Create a minimal report structure
                                backup_data = {
                                    "title": "PharmCAT Report (Backup)",
                                    "timestamp": datetime.datetime.now().isoformat(),
                                    "pharmcatVersion": "2.15.5",
                                    "dataVersion": "N/A",
                                    "genes": phenotype_data.get("phenotypes", {}),
                                    "drugs": {},
                                    "messages": [{"text": "This is a backup report created for inspection"}]
                                }
                            else:
                                # Create a minimal empty structure
                                backup_data = {
                                    "title": "PharmCAT Report (Backup)",
                                    "timestamp": datetime.datetime.now().isoformat(),
                                    "pharmcatVersion": "2.15.5",
                                    "dataVersion": "N/A",
                                    "genes": {},
                                    "drugs": {},
                                    "messages": [{"text": "This is a backup report created for inspection. Original report.json was not found."}]
                                }
                                
                            # Write the backup JSON file
                            with open(backup_json_path, 'w') as f:
                                json.dump(backup_data, f, indent=2)
                            logger.info(f"Backup JSON report created at {backup_json_path}")
                        except Exception as e:
                            logger.error(f"Error creating backup JSON report: {str(e)}")
                
                # Process report.json file - our single source of truth
                report_json_path = Path(temp_dir) / report_json_file
                
                # Read and parse the report JSON file
                with open(report_json_path, 'r') as f:
                    report_data = json.load(f)
                
                logger.info(f"Loaded report data successfully. Keys: {list(report_data.keys())}")
                
                # Always create a permanent copy of the raw report.json in the patient directory
                # Try a more direct approach to ensure the file is created
                raw_report_path = patient_dir / f"{base_name}_raw_report.json"
                try:
                    with open(raw_report_path, 'w') as f:
                        json.dump(report_data, f, indent=2)
                    logger.info(f"Raw report.json saved to {raw_report_path}")
                    
                    # Double check the file exists
                    if os.path.exists(raw_report_path):
                        logger.info(f"Verified raw report.json file exists at {raw_report_path} with size {os.path.getsize(raw_report_path)}")
                    else:
                        logger.error(f"Failed to create raw report.json at {raw_report_path}")
                except Exception as e:
                    logger.error(f"Error saving raw report.json: {str(e)}")
                    
                # We no longer create backup copies since patient directories work better
                # backup_report_path = reports_dir / "latest_pharmcat_report.json"
                # try:
                #     with open(backup_report_path, 'w') as f:
                #         json.dump(report_data, f, indent=2)
                #     logger.info(f"Backup report.json saved to {backup_report_path}")
                
                # Create a standard JSON report in the patient directory
                try:
                    standard_json_path = patient_dir / f"{base_name}_pgx_report.json"
                    with open(standard_json_path, 'w') as f:
                        json.dump(report_data, f, indent=2)
                    logger.info(f"Standard JSON report saved to {standard_json_path}")
                except Exception as e:
                    logger.error(f"Error saving standard JSON report: {str(e)}")
                
                # Return success response with report URLs
                return jsonify({
                    "success": True,
                    "message": "PharmCAT analysis completed successfully",
                    "data": {
                        "job_id": base_name,
                        "html_report_url": f"/reports/{base_name}/{base_name}_pgx_report.html" if report_html_file else None,
                        "json_report_url": f"/reports/{base_name}/{base_name}_pgx_report.json" if report_json_file else None,
                        "tsv_report_url": f"/reports/{base_name}/{base_name}_pgx_report.tsv" if report_tsv_file else None,
                        "raw_report_url": f"/reports/{base_name}/{base_name}_raw_report.json",
                        # Add standard URLs for our custom reports
                        "pdf_report_url": f"/reports/{base_name}/{base_name}_pgx_report.pdf",
                        "interactive_html_report_url": f"/reports/{base_name}/{base_name}_pgx_report_interactive.html",
                        # Add normalized URLs for PharmCAT original reports
                        "pharmcat_html_report_url": f"/reports/{base_name}/{base_name}_pgx_pharmcat.html" if report_html_file else None,
                        "pharmcat_json_report_url": f"/reports/{base_name}/{base_name}_pgx_pharmcat.json" if report_json_file else None,
                        "pharmcat_tsv_report_url": f"/reports/{base_name}/{base_name}_pgx_pharmcat.tsv" if report_tsv_file else None,
                        "genes": report_data.get("genes", []),
                        "drugs": report_data.get("drugs", []),
                        "messages": report_data.get("messages", [])
                    }
                })
            
            except subprocess.CalledProcessError as e:
                logger.error(f"PharmCAT process error: {e.stderr}")
                return jsonify({
                    "success": False,
                    "message": f"PharmCAT process error: {e.stderr}"
                }), 500
                
            except Exception as e:
                logger.error(f"Error running PharmCAT: {str(e)}")
                logger.error(traceback.format_exc())
                return jsonify({
                    "success": False,
                    "message": f"Error running PharmCAT: {str(e)}"
                }), 500
                
    except Exception as e:
        logger.error(f"Error processing request: {str(e)}")
        logger.error(traceback.format_exc())
        return jsonify({
            "success": False,
            "message": f"Error: {str(e)}"
        }), 500

@app.route('/test', methods=['GET'])
def test_connection():
    """
    Test endpoint to verify connectivity.
    """
    return jsonify({
        "status": "ok",
        "message": "PharmCAT wrapper is running",
        "data_dir": os.getenv('DATA_DIR', '/data'),
        "pharmcat_jar": PHARMCAT_JAR,
        "jar_exists": os.path.exists(PHARMCAT_JAR)
    })

def run_pharmcat_jar(input_file: str, output_dir: str, sample_id: Optional[str] = None, report_id: Optional[str] = None) -> Dict[str, Any]:
    """
    Run PharmCAT directly using the pipeline distribution.
    
    Args:
        input_file: Path to the input file
        output_dir: Directory to store the results
        sample_id: Optional sample ID to use
        report_id: Optional report ID to use for consistent directory naming
        
    Returns:
        Dictionary containing PharmCAT results
    """
    try:
        # Make sure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        # Use report_id if provided, otherwise use filename without extension as base name
        if report_id:
            base_name = report_id
            logger.info(f"Using provided report ID: {base_name}")
        else:
            # Use just the filename without extension as base name
            base_name = Path(input_file).stem
            if base_name.endswith('.vcf'):
                base_name = base_name[:-4]
            logger.info(f"Using file stem as base name: {base_name}")
        
        # Use pharmcat_pipeline as per official documentation
        logger.info(f"Running PharmCAT pipeline with input: {input_file}")
        
        # Prepare command
        cmd = [
            "pharmcat_pipeline",
            "-G",  # Bypass gVCF check
            "-o", output_dir,  # Output directory
            "-v",  # Verbose output
            "-reporterHtml",    # Generate HTML report
            "-reporterJson",    # Generate JSON report 
            "-reporterCallsOnlyTsv",  # Generate TSV report with only calls data
            input_file  # Input file should be the last argument
        ]
        
        # Add sample ID if provided
        if sample_id:
            cmd.extend(["-s", sample_id])
        
        # Set environment variables
        env = os.environ.copy()
        env["JAVA_TOOL_OPTIONS"] = "-Xmx4g -XX:+UseG1GC"
        env["PHARMCAT_LOG_LEVEL"] = "DEBUG"
        
        # Log important environment info for debugging
        logger.info(f"PHARMCAT_JAR location: {PHARMCAT_JAR}")
        logger.info(f"PHARMCAT_PIPELINE_DIR: {PHARMCAT_PIPELINE_DIR}")
        logger.info(f"PATH environment: {env.get('PATH', 'Not set')}")
        
        # Check if pharmcat_pipeline exists and is executable
        try:
            subprocess.run(["which", "pharmcat_pipeline"], check=True, capture_output=True, text=True)
            logger.info("pharmcat_pipeline command found in PATH")
        except subprocess.CalledProcessError:
            logger.warning("pharmcat_pipeline command NOT found in PATH")
        
        # Run PharmCAT pipeline
        logger.info(f"Executing PharmCAT command: {' '.join(cmd)}")
        process = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
            env=env,
            timeout=300  # 5 minute timeout
        )
        
        logger.info(f"PharmCAT execution completed with output: {process.stdout}")
        
        # Expected output files
        results_files = {
            "match_results": os.path.join(output_dir, f"{base_name}.match.json"),
            "phenotype_results": os.path.join(output_dir, f"{base_name}.phenotype.json"),
            "report_html": os.path.join(output_dir, f"{base_name}.report.html"),
            "report_json": os.path.join(output_dir, f"{base_name}.report.json"),
            "report_tsv": os.path.join(output_dir, f"{base_name}.report.tsv")
        }
        
        # Check if essential files exist
        essential_files = ["match_results", "phenotype_results"]
        if not all(os.path.exists(results_files[f]) for f in essential_files):
            missing_files = [k for k in essential_files if not os.path.exists(results_files[k])]
            error_msg = f"Missing required output files: {', '.join(missing_files)}"
            logger.error(error_msg)
            return {
                "success": False,
                "message": error_msg
            }
        
        # Read results
        results = {}
        for name, path in results_files.items():
            if os.path.exists(path):
                with open(path, 'r') as f:
                    if name in ["report_html", "report_tsv"]:
                        results[name] = f.read()
                    else:
                        results[name] = json.load(f)
            else:
                logger.warning(f"Optional output file not found: {path}")
        
        # Copy PharmCAT reports to /data/reports for direct access
        reports_dir = Path("/data/reports")
        reports_dir.mkdir(parents=True, exist_ok=True)
        
        # Create a patient-specific directory for this job
        patient_dir = reports_dir / base_name
        patient_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Created patient-specific directory: {patient_dir}")
        
        # Map of source files to destination files
        report_files = {
            f"{base_name}.report.html": f"{base_name}_pgx_report.html",
            f"{base_name}.report.json": f"{base_name}_pgx_report.json",
            f"{base_name}.report.tsv": f"{base_name}_pgx_report.tsv",
            f"{base_name}.match.json": f"{base_name}_pgx_match.json",
            f"{base_name}.phenotype.json": f"{base_name}_pgx_phenotype.json"
        }
        
        # Copy all report files that exist
        for src_name, dest_name in report_files.items():
            src_path = Path(output_dir) / src_name
            dest_path = patient_dir / dest_name
            
            if os.path.exists(src_path):
                shutil.copy2(src_path, dest_path)
                logger.info(f"Report file copied to {dest_path}")
            else:
                logger.warning(f"Report file not found at {src_path}")
                
        # Extract gene data and drug recommendations
        genes_data = []
        drug_recommendations = []
        
        phenotype_path = Path(output_dir) / f"{base_name}.phenotype.json"
        if os.path.exists(phenotype_path):
            try:
                with open(phenotype_path, 'r') as f:
                    phenotype_data = json.load(f)
                    
                # Extract gene data from phenotype file
                if "phenotypes" in phenotype_data:
                    for gene_id, gene_info in phenotype_data["phenotypes"].items():
                        gene_entry = {
                            "gene": gene_id,
                            "diplotype": {
                                "name": gene_info.get("diplotype", "Unknown"),
                                "activityScore": gene_info.get("activityScore")
                            },
                            "phenotype": {
                                "info": gene_info.get("phenotype", "Unknown")
                            }
                        }
                        genes_data.append(gene_entry)
                        
                # Extract drug recommendations
                if "drugRecommendations" in phenotype_data:
                    drug_recommendations = phenotype_data["drugRecommendations"]
                    
                logger.info(f"Extracted {len(genes_data)} genes and {len(drug_recommendations)} drug recommendations")
            except Exception as e:
                logger.error(f"Error parsing phenotype file: {str(e)}")
        
        # Prepare the result data
        return {
            "success": True,
            "message": "PharmCAT analysis completed successfully",
            "data": {
                "job_id": base_name,
                "pdf_report_url": f"/reports/{base_name}/{base_name}_pgx_report.pdf",
                "html_report_url": f"/reports/{base_name}/{base_name}_pgx_report.html",
                "interactive_html_report_url": f"/reports/{base_name}/{base_name}_pgx_report_interactive.html",
                "json_report_url": f"/reports/{base_name}/{base_name}_pgx_report.json",
                "tsv_report_url": f"/reports/{base_name}/{base_name}_pgx_report.tsv",
                "match_json_url": f"/reports/{base_name}/{base_name}_pgx_match.json",
                "phenotype_json_url": f"/reports/{base_name}/{base_name}_pgx_phenotype.json",
                # Add normalized URLs for PharmCAT original reports
                "pharmcat_html_report_url": f"/reports/{base_name}/{base_name}_pgx_pharmcat.html",
                "pharmcat_json_report_url": f"/reports/{base_name}/{base_name}_pgx_pharmcat.json",
                "pharmcat_tsv_report_url": f"/reports/{base_name}/{base_name}_pgx_pharmcat.tsv",
                "genes": genes_data,
                "drugRecommendations": drug_recommendations,
                "results": results
            }
        }
    
    except subprocess.CalledProcessError as e:
        error_msg = f"PharmCAT execution failed: {e.stderr}" if e.stderr else str(e)
        logger.error(error_msg)
        return {
            "success": False,
            "message": "PharmCAT execution failed",
            "error": str(e),
            "stderr": e.stderr if hasattr(e, 'stderr') else None
        }
    
    except Exception as e:
        error_msg = f"Error running PharmCAT: {str(e)}"
        logger.error(error_msg)
        return {
            "success": False,
            "message": error_msg
        }

if __name__ == '__main__':
    # Start the Flask app
    logger.info("Starting PharmCAT wrapper service with DATA_DIR=%s", os.getenv('DATA_DIR', '/data'))
    logger.info("PharmCAT JAR location: %s", PHARMCAT_JAR)
    
    # Configure Flask to listen on all interfaces
    app.run(
        host='0.0.0.0',
        port=5000,
        debug=True,
        use_reloader=False  # Disable reloader in debug mode to avoid duplicate processes
    ) 