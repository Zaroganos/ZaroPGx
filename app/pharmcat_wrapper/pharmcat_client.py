import os
import requests
import logging
import subprocess
import json
import tempfile
from typing import Dict, List, Any, Optional

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# PharmCAT service configuration
PHARMCAT_SERVICE_URL = os.environ.get("PHARMCAT_SERVICE_URL", "http://pharmcat:8080")
PHARMCAT_DOCKER_IMAGE = os.environ.get("PHARMCAT_DOCKER_IMAGE", "pgkb/pharmcat:latest")
PHARMCAT_JAR_PATH = os.environ.get("PHARMCAT_JAR_PATH", "/app/lib/pharmcat.jar")

def call_pharmcat_service(input_file: str) -> Dict[str, Any]:
    """
    Call the PharmCAT service to process a genetic data file.
    
    Args:
        input_file: Path to the input file (23andMe format)
        
    Returns:
        Dictionary containing PharmCAT results
    """
    try:
        logger.info(f"Calling PharmCAT service for file: {input_file}")
        
        # Try direct JAR execution or wrapper API
        pharmcat_jar = os.environ.get("PHARMCAT_JAR_PATH", "/pharmcat/pharmcat.jar")
        
        # First check if we can access the JAR directly
        if os.path.exists(pharmcat_jar):
            logger.info(f"Found PharmCAT JAR at {pharmcat_jar}, using direct execution")
            results = run_pharmcat_jar(input_file)
        else:
            # Try the wrapper API instead
            logger.info("PharmCAT JAR not found, trying wrapper API")
            pharmcat_api_url = os.environ.get("PHARMCAT_API_URL", "http://pharmcat-wrapper:5000")
            with open(input_file, 'rb') as f:
                files = {'file': f}
                response = requests.post(
                    f"{pharmcat_api_url}/process",
                    files=files,
                    timeout=300  # 5 minute timeout
                )
                response.raise_for_status()
                results = response.json()
            
        return results
    except Exception as e:
        logger.error(f"Error calling PharmCAT: {str(e)}")
        raise

def call_pharmcat_api(input_file: str) -> Dict[str, Any]:
    """
    Call PharmCAT REST API service.
    
    Args:
        input_file: Path to the input file
        
    Returns:
        Dictionary containing PharmCAT results
    """
    try:
        # Open the file for reading
        with open(input_file, 'rb') as f:
            # Prepare the file for upload
            files = {'file': f}
            
            # Make the POST request
            response = requests.post(
                f"{PHARMCAT_SERVICE_URL}/api/process",
                files=files,
                timeout=300  # 5 minute timeout
            )
            
            # Check if request was successful
            response.raise_for_status()
            
            # Parse response
            results = response.json()
            logger.info(f"PharmCAT API call successful: {len(results)} results returned")
            
            return results
    except requests.exceptions.RequestException as e:
        logger.error(f"Error calling PharmCAT API: {str(e)}")
        raise

def run_pharmcat_jar(input_file: str) -> Dict[str, Any]:
    """
    Run PharmCAT directly using the JAR file.
    
    Args:
        input_file: Path to the input file
        
    Returns:
        Dictionary containing PharmCAT results
    """
    try:
        # Create a temporary directory for output
        with tempfile.TemporaryDirectory() as temp_dir:
            # Get the JAR path from environment or use default
            pharmcat_jar = os.environ.get("PHARMCAT_JAR_PATH", "/pharmcat/pharmcat.jar")
            logger.info(f"Running PharmCAT JAR on file: {input_file}, output dir: {temp_dir}")
            
            # Execute PharmCAT JAR
            command = [
                "java", "-jar", pharmcat_jar,
                "-m", "vcf",  # Use VCF mode instead of 23andme
                "-i", input_file,
                "-o", temp_dir
            ]
            
            process = subprocess.run(
                command,
                capture_output=True,
                text=True,
                check=True
            )
            
            logger.info(f"PharmCAT execution successful: {process.stdout}")
            
            # Read the results JSON file
            results_file = os.path.join(temp_dir, "results.json")
            if os.path.exists(results_file):
                with open(results_file, 'r') as f:
                    results = json.load(f)
                return results
            else:
                raise FileNotFoundError(f"PharmCAT results file not found: {results_file}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running PharmCAT: {e.stderr}")
        raise
    except Exception as e:
        logger.error(f"Error processing PharmCAT results: {str(e)}")
        raise

def run_pharmcat_docker(input_file: str) -> Dict[str, Any]:
    """
    Run PharmCAT using Docker container.
    
    Args:
        input_file: Path to the input file
        
    Returns:
        Dictionary containing PharmCAT results
    """
    try:
        # Create a temporary directory for output
        with tempfile.TemporaryDirectory() as temp_dir:
            # Get the absolute path to the input file and temporary directory
            abs_input_file = os.path.abspath(input_file)
            abs_temp_dir = os.path.abspath(temp_dir)
            
            # Extract directory and filename
            input_dir = os.path.dirname(abs_input_file)
            input_filename = os.path.basename(abs_input_file)
            
            logger.info(f"Running PharmCAT Docker for file: {input_filename}")
            
            # Docker run command
            command = [
                "docker", "run", "--rm",
                "-v", f"{input_dir}:/input",
                "-v", f"{abs_temp_dir}:/output",
                PHARMCAT_DOCKER_IMAGE,
                "-m", "23andme",
                "-i", f"/input/{input_filename}",
                "-o", "/output"
            ]
            
            process = subprocess.run(
                command,
                capture_output=True,
                text=True,
                check=True
            )
            
            logger.info(f"PharmCAT Docker execution successful: {process.stdout}")
            
            # Read the results JSON file
            results_file = os.path.join(temp_dir, "results.json")
            if os.path.exists(results_file):
                with open(results_file, 'r') as f:
                    results = json.load(f)
                return results
            else:
                raise FileNotFoundError(f"PharmCAT results file not found: {results_file}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running PharmCAT Docker: {e.stderr}")
        raise
    except Exception as e:
        logger.error(f"Error processing PharmCAT Docker results: {str(e)}")
        raise

def parse_pharmcat_results(results: Dict[str, Any]) -> Dict[str, Any]:
    """
    Parse PharmCAT results into a standardized format.
    
    Args:
        results: Raw PharmCAT results
        
    Returns:
        Standardized results dictionary
    """
    try:
        # Extract gene-specific information
        genes = {}
        for gene in results.get("genes", []):
            gene_symbol = gene.get("gene")
            diplotype = gene.get("diplotype", {}).get("name", "Unknown")
            phenotype = gene.get("phenotype", {}).get("info", "Unknown")
            activity_score = gene.get("diplotype", {}).get("activityScore")
            
            genes[gene_symbol] = {
                "diplotype": diplotype,
                "phenotype": phenotype,
                "activity_score": activity_score
            }
        
        # Extract drug recommendations
        recommendations = []
        for drug in results.get("drugRecommendations", []):
            drug_name = drug.get("drug", {}).get("name", "Unknown")
            recommendations.append({
                "drug": drug_name,
                "guidelines": drug.get("guidelines", []),
                "recommendation": drug.get("recommendation", "Unknown"),
                "classification": drug.get("classification", "Unknown")
            })
        
        # Return standardized format
        return {
            "genes": genes,
            "recommendations": recommendations,
            "report_id": results.get("reportId"),
            "report_time": results.get("reportTime")
        }
    except Exception as e:
        logger.error(f"Error parsing PharmCAT results: {str(e)}")
        raise 