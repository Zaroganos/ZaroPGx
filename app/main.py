from fastapi import FastAPI, Depends, HTTPException, status, Request, File, UploadFile, Form
from fastapi.middleware.cors import CORSMiddleware
from fastapi.security import OAuth2PasswordBearer, OAuth2PasswordRequestForm
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from fastapi.responses import HTMLResponse, RedirectResponse
from jose import JWTError, jwt
from datetime import datetime, timedelta
from typing import Optional
import os
import shutil
from pathlib import Path
from dotenv import load_dotenv
import logging

from app.api.routes import upload_router, report_router
from app.api.models import Token, TokenData

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Security configuration
SECRET_KEY = os.getenv("SECRET_KEY", "supersecretkey")  # In production, use env var
ALGORITHM = "HS256"
ACCESS_TOKEN_EXPIRE_MINUTES = 30

# Directory setup
BASE_DIR = Path(__file__).resolve().parent
TEMPLATE_DIR = BASE_DIR / "templates"
UPLOAD_DIR = Path("/data/uploads")
UPLOAD_DIR.mkdir(parents=True, exist_ok=True)

# Initialize templates
templates = Jinja2Templates(directory=str(TEMPLATE_DIR))

# Initialize FastAPI app
app = FastAPI(
    title="ZaroPGx - Pharmacogenomic Report Generator",
    description="API for processing genetic data and generating pharmacogenomic reports",
    version="0.1.0"
)

# OAuth2
oauth2_scheme = OAuth2PasswordBearer(tokenUrl="token")

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # In production, specify exact origins
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Include routers
app.include_router(upload_router.router)
app.include_router(report_router.router)

# JWT token functions
def create_access_token(data: dict, expires_delta: Optional[timedelta] = None):
    to_encode = data.copy()
    if expires_delta:
        expire = datetime.utcnow() + expires_delta
    else:
        expire = datetime.utcnow() + timedelta(minutes=15)
    to_encode.update({"exp": expire})
    encoded_jwt = jwt.encode(to_encode, SECRET_KEY, algorithm=ALGORITHM)
    return encoded_jwt

async def get_current_user(token: str = Depends(oauth2_scheme)):
    credentials_exception = HTTPException(
        status_code=status.HTTP_401_UNAUTHORIZED,
        detail="Could not validate credentials",
        headers={"WWW-Authenticate": "Bearer"},
    )
    try:
        payload = jwt.decode(token, SECRET_KEY, algorithms=[ALGORITHM])
        username: str = payload.get("sub")
        if username is None:
            raise credentials_exception
        token_data = TokenData(username=username)
    except JWTError:
        raise credentials_exception
    # In a real app, get user from database
    if token_data.username != "testuser":  # Mock user validation
        raise credentials_exception
    return token_data.username

# Authentication endpoint
@app.post("/token", response_model=Token)
async def login_for_access_token(form_data: OAuth2PasswordRequestForm = Depends()):
    # In a real app, validate against database
    if form_data.username != "testuser" or form_data.password != "testpassword":
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Incorrect username or password",
            headers={"WWW-Authenticate": "Bearer"},
        )
    access_token_expires = timedelta(minutes=ACCESS_TOKEN_EXPIRE_MINUTES)
    access_token = create_access_token(
        data={"sub": form_data.username}, expires_delta=access_token_expires
    )
    return {"access_token": access_token, "token_type": "bearer"}

@app.get("/", response_class=HTMLResponse)
async def home(request: Request):
    return templates.TemplateResponse("index.html", {"request": request})

@app.get("/api")
async def api_root():
    return {"message": "Welcome to ZaroPGx API", "docs": "/docs"}

# Make the health check endpoint simple and dependency-free
@app.get("/health")
async def health_check():
    logger.info("Health check called")
    return {"status": "healthy", "timestamp": str(datetime.utcnow())}

@app.post("/upload-vcf")
async def upload_vcf_file(
    request: Request,
    vcfFile: UploadFile = File(...),
    sampleId: str = Form(None)
):
    logger.info(f"Received VCF file upload: {vcfFile.filename}, Sample ID: {sampleId}")
    
    # Generate a unique filename if no sample ID provided
    filename = f"{sampleId or 'sample'}-{datetime.now().strftime('%Y%m%d%H%M%S')}.vcf"
    file_path = UPLOAD_DIR / filename
    
    # Save the uploaded file
    try:
        with open(file_path, "wb") as buffer:
            shutil.copyfileobj(vcfFile.file, buffer)
        
        # In a real implementation, you would call the PharmCAT and Aldy services here
        # For now, just return a success message
        return templates.TemplateResponse(
            "index.html", 
            {
                "request": request, 
                "message": f"File {vcfFile.filename} uploaded successfully! Processing will begin shortly.",
                "success": True
            }
        )
    except Exception as e:
        logger.error(f"Error processing file upload: {str(e)}")
        return templates.TemplateResponse(
            "index.html", 
            {
                "request": request, 
                "message": f"Error processing file: {str(e)}",
                "success": False
            }
        ) 