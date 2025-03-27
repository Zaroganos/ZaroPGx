import os
import json
import logging
from datetime import datetime
from typing import List, Dict, Any
from weasyprint import HTML, CSS
from jinja2 import Environment, FileSystemLoader

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Constants
TEMPLATE_DIR = os.path.join(os.path.dirname(__file__), "templates")
CSS_FILE = os.path.join(TEMPLATE_DIR, "style.css")

# Initialize Jinja2 environment
env = Environment(loader=FileSystemLoader(TEMPLATE_DIR))

def generate_pdf_report(
    patient_id: str,
    report_id: str,
    diplotypes: List[Dict[str, Any]],
    recommendations: List[Dict[str, Any]],
    report_path: str
) -> str:
    """
    Generate a PDF pharmacogenomic report.
    
    Args:
        patient_id: Patient identifier
        report_id: Report identifier
        diplotypes: List of diplotype results
        recommendations: List of drug recommendations
        report_path: Path to save the PDF report
        
    Returns:
        Path to the generated PDF file
    """
    try:
        logger.info(f"Generating PDF report for patient {patient_id}")
        
        # Ensure the directory exists
        os.makedirs(os.path.dirname(report_path), exist_ok=True)
        
        # Prepare the report data
        report_data = {
            "patient_id": patient_id,
            "report_id": report_id,
            "report_date": datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC"),
            "diplotypes": diplotypes,
            "recommendations": recommendations,
            "disclaimer": get_disclaimer()
        }
        
        # Load and render the HTML template
        template = env.get_template("report_template.html")
        html_content = template.render(**report_data)
        
        # Generate PDF from HTML
        generate_pdf_from_html(html_content, report_path)
        
        logger.info(f"PDF report generated successfully: {report_path}")
        return report_path
    except Exception as e:
        logger.error(f"Error generating PDF report: {str(e)}")
        raise

def generate_pdf_from_html(html_content: str, output_path: str) -> None:
    """
    Generate PDF from HTML content using WeasyPrint.
    
    Args:
        html_content: HTML content to convert
        output_path: Path to save the PDF
    """
    try:
        # Load CSS if it exists
        css = None
        if os.path.exists(CSS_FILE):
            css = CSS(filename=CSS_FILE)
        
        # Generate PDF
        html = HTML(string=html_content)
        if css:
            html.write_pdf(output_path, stylesheets=[css])
        else:
            html.write_pdf(output_path)
    except Exception as e:
        logger.error(f"Error converting HTML to PDF: {str(e)}")
        raise

def get_disclaimer() -> str:
    """
    Return the legal disclaimer for pharmacogenomic reports.
    """
    return """
    DISCLAIMER: This pharmacogenomic report is for informational purposes only and is not intended
    to be used as a substitute for professional medical advice, diagnosis, or treatment. The recommendations
    in this report are based on guidelines from the Clinical Pharmacogenetics Implementation Consortium (CPIC)
    and are subject to change as new research becomes available.
    
    The results should be interpreted by a healthcare professional in the context of the patient's
    clinical situation. Medication decisions should never be made solely based on this report without
    consulting a qualified healthcare provider.
    """

def organize_gene_drug_recommendations(recommendations: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    """
    Organize recommendations by gene and drug for better presentation.
    
    Args:
        recommendations: List of drug recommendations
        
    Returns:
        Dictionary organized by gene and drug
    """
    organized = {}
    
    for rec in recommendations:
        gene = rec.get("gene")
        drug = rec.get("drug")
        
        if gene not in organized:
            organized[gene] = {}
            
        if drug not in organized[gene]:
            organized[gene][drug] = []
            
        organized[gene][drug].append(rec)
    
    return organized

def create_interactive_html_report(
    patient_id: str,
    report_id: str,
    diplotypes: List[Dict[str, Any]],
    recommendations: List[Dict[str, Any]],
    output_path: str
) -> str:
    """
    Create an interactive HTML report with JavaScript visualizations.
    
    Args:
        patient_id: Patient identifier
        report_id: Report identifier
        diplotypes: List of diplotype results
        recommendations: List of drug recommendations
        output_path: Path to save the HTML report
        
    Returns:
        Path to the generated HTML file
    """
    try:
        logger.info(f"Generating interactive HTML report for patient {patient_id}")
        
        # Ensure the directory exists
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        # Prepare the report data
        report_data = {
            "patient_id": patient_id,
            "report_id": report_id,
            "report_date": datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC"),
            "diplotypes": diplotypes,
            "recommendations": recommendations,
            "organized_recommendations": json.dumps(organize_gene_drug_recommendations(recommendations)),
            "disclaimer": get_disclaimer()
        }
        
        # Load and render the HTML template
        template = env.get_template("interactive_report.html")
        html_content = template.render(**report_data)
        
        # Write HTML to file
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(html_content)
        
        logger.info(f"Interactive HTML report generated successfully: {output_path}")
        return output_path
    except Exception as e:
        logger.error(f"Error generating interactive HTML report: {str(e)}")
        raise 