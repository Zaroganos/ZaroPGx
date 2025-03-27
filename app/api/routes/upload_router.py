import os
import uuid
import shutil
from fastapi import APIRouter, Depends, UploadFile, File, Form, HTTPException, BackgroundTasks
from sqlalchemy.orm import Session
from typing import Optional
import logging
from datetime import datetime

from app.api.db import get_db, create_patient, register_genetic_data
from app.api.models import UploadResponse, FileType
from app.pharmcat_wrapper.pharmcat_client import call_pharmcat_service
from ..utils.security import get_current_user

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Initialize router
router = APIRouter(
    prefix="/upload",
    tags=["upload"],
    dependencies=[Depends(get_current_user)]
)

# Constants
UPLOAD_DIR = os.environ.get("UPLOAD_DIR", "/app/data/uploads")
# Ensure upload directory exists
os.makedirs(UPLOAD_DIR, exist_ok=True)

# Background task to process 23andMe data
def process_23andme_background(file_path: str, patient_id: str, data_id: str):
    try:
        logger.info(f"Processing 23andMe data for patient {patient_id}, file {data_id}")
        # Call PharmCAT service for processing
        results = call_pharmcat_service(file_path)
        # Update database with results
        logger.info(f"23andMe processing complete: {results}")
    except Exception as e:
        logger.error(f"Error processing 23andMe data: {str(e)}")
        # Update status in database to failed
        # This would be implemented with a db call

# Background task to process VCF data
def process_vcf_background(file_path: str, patient_id: str, data_id: str):
    try:
        logger.info(f"Processing VCF data for patient {patient_id}, file {data_id}")
        # Call Aldy service for processing
        # This would be implemented with a call to Aldy
        logger.info(f"VCF processing complete")
    except Exception as e:
        logger.error(f"Error processing VCF data: {str(e)}")
        # Update status in database to failed

# File upload endpoints
@router.post("/23andme", response_model=UploadResponse)
async def upload_23andme(
    background_tasks: BackgroundTasks,
    file: UploadFile = File(...),
    patient_identifier: Optional[str] = Form(None),
    db: Session = Depends(get_db)
):
    """
    Upload 23andMe genetic data file for processing.
    Returns a unique ID for tracking the processing status.
    """
    # Validate file
    if not file.filename.endswith(('.txt', '.csv')):
        raise HTTPException(status_code=400, detail="Invalid file format. 23andMe files should be text files.")
    
    # Generate unique identifiers
    file_id = str(uuid.uuid4())
    if not patient_identifier:
        patient_identifier = f"patient_{uuid.uuid4()}"
    
    try:
        # Create patient record
        patient_id = create_patient(db, patient_identifier)
        
        # Create directory for patient if it doesn't exist
        patient_dir = os.path.join(UPLOAD_DIR, str(patient_id))
        os.makedirs(patient_dir, exist_ok=True)
        
        # Save file with secure name
        file_path = os.path.join(patient_dir, f"{file_id}.txt")
        with open(file_path, "wb") as buffer:
            shutil.copyfileobj(file.file, buffer)
        
        # Register genetic data in database
        data_id = register_genetic_data(db, patient_id, FileType.TWENTYTHREE_AND_ME.value, file_path)
        
        # Schedule background processing
        background_tasks.add_task(process_23andme_background, file_path, str(patient_id), str(data_id))
        
        return UploadResponse(
            file_id=str(data_id),
            file_type=FileType.TWENTYTHREE_AND_ME,
            status="queued",
            message="File uploaded successfully. Processing started."
        )
        
    except Exception as e:
        logger.error(f"Error during 23andMe upload: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Error processing upload: {str(e)}")

@router.post("/wgs", response_model=UploadResponse)
async def upload_wgs(
    background_tasks: BackgroundTasks,
    file: UploadFile = File(...),
    patient_identifier: Optional[str] = Form(None),
    db: Session = Depends(get_db)
):
    """
    Upload whole genome sequencing (WGS) VCF file for processing.
    Returns a unique ID for tracking the processing status.
    """
    # Validate file
    if not file.filename.endswith('.vcf'):
        raise HTTPException(status_code=400, detail="Invalid file format. WGS files should be in VCF format.")
    
    # Generate unique identifiers
    file_id = str(uuid.uuid4())
    if not patient_identifier:
        patient_identifier = f"patient_{uuid.uuid4()}"
    
    try:
        # Create patient record
        patient_id = create_patient(db, patient_identifier)
        
        # Create directory for patient if it doesn't exist
        patient_dir = os.path.join(UPLOAD_DIR, str(patient_id))
        os.makedirs(patient_dir, exist_ok=True)
        
        # Save file with secure name
        file_path = os.path.join(patient_dir, f"{file_id}.vcf")
        with open(file_path, "wb") as buffer:
            shutil.copyfileobj(file.file, buffer)
        
        # Register genetic data in database
        data_id = register_genetic_data(db, patient_id, FileType.VCF.value, file_path)
        
        # Schedule background processing
        background_tasks.add_task(process_vcf_background, file_path, str(patient_id), str(data_id))
        
        return UploadResponse(
            file_id=str(data_id),
            file_type=FileType.VCF,
            status="queued",
            message="File uploaded successfully. Processing started."
        )
        
    except Exception as e:
        logger.error(f"Error during VCF upload: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Error processing upload: {str(e)}")

@router.get("/status/{file_id}")
async def get_upload_status(file_id: str, db: Session = Depends(get_db)):
    """
    Check the processing status of an uploaded file.
    """
    try:
        # Query processing status from database
        # This would be implemented with a db query
        # For now, return mock data
        return {
            "file_id": file_id,
            "status": "processing",
            "progress": 50,
            "message": "Analysis in progress"
        }
    except Exception as e:
        logger.error(f"Error getting status: {str(e)}")
        raise HTTPException(status_code=404, detail="File not found") 