from pydantic import BaseModel, Field
from typing import List, Dict, Optional, Any
from datetime import datetime
from enum import Enum


class Token(BaseModel):
    access_token: str
    token_type: str


class TokenData(BaseModel):
    username: Optional[str] = None


class FileType(str, Enum):
    TWENTYTHREE_AND_ME = "23andme"
    VCF = "vcf"
    OTHER = "other"


class UploadResponse(BaseModel):
    file_id: str
    file_type: FileType
    status: str
    message: str


class ProcessingStatus(str, Enum):
    QUEUED = "queued"
    PROCESSING = "processing"
    COMPLETED = "completed"
    FAILED = "failed"


class GeneticDataStatus(BaseModel):
    file_id: str
    file_type: FileType
    status: ProcessingStatus
    created_at: datetime
    processed_at: Optional[datetime] = None
    error_message: Optional[str] = None


class Allele(BaseModel):
    name: str
    function: Optional[str] = None
    activity_score: Optional[float] = None


class Diplotype(BaseModel):
    gene: str
    diplotype: str
    phenotype: Optional[str] = None
    activity_score: Optional[float] = None
    confidence: Optional[float] = None
    calling_method: str


class AlleleCallResult(BaseModel):
    patient_id: str
    file_id: str
    diplotypes: List[Diplotype]
    created_at: datetime


class DrugRecommendation(BaseModel):
    drug: str
    gene: str
    guideline: str
    recommendation: str
    classification: str  # e.g., "Strong", "Moderate"
    literature_references: Optional[List[str]] = None


class ReportRequest(BaseModel):
    patient_id: str
    file_id: str
    report_type: str = "comprehensive"
    include_drugs: Optional[List[str]] = None  # If None, include all drugs


class ReportResponse(BaseModel):
    report_id: str
    patient_id: str
    created_at: datetime
    report_url: str
    report_type: str 