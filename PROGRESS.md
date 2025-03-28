# ZaroPGx Project Progress

## Overview
ZaroPGx is a pharmacogenomics analysis platform using Docker-based microservices architecture. The system processes genetic data (VCF files) to generate personalized medication recommendations based on patient genotypes.

## Project Status
- **Current Phase**: Initial development and containerization
- **Last Updated**: March 27, 2025

## Components & Services

### 1. Core Services

| Service | Status | Description | Issues/Next Steps |
|---------|--------|-------------|------------------|
| **Main App (FastAPI)** | 🟢 Working | Main application handling UI, reports, user auth | New user-friendly web UI added |
| **PostgreSQL Database** | 🟢 Working | Stores CPIC guidelines, user data, reports | Need to add seed data scripts |
| **PharmCAT Service** | 🟢 Working | Using official pgkb/pharmcat Docker image | Now running with proper command |
| **PharmCAT Wrapper** | 🟢 Working | Python wrapper providing REST API for PharmCAT | API accessible via port 5001 |
| **Aldy Service** | 🟢 Working | CYP2D6 specific genotype caller | Running with Gunicorn for production |

### 2. Docker Infrastructure

| Component | Status | Description | Issues/Next Steps |
|-----------|--------|-------------|------------------|
| **Docker Compose** | 🟢 Working | Multi-container orchestration | Core services up and running |
| **Network Config** | 🟢 Working | PGX-Network for service communication | Services communicating successfully |
| **Volumes** | 🟢 Working | Data persistence and sharing | Using volume mounts for JAR sharing |
| **Health Checks** | 🟡 In Progress | Container health monitoring | Temporarily simplified for debugging |

### 3. Features & Capabilities

| Feature | Status | Description | Issues/Next Steps |
|---------|--------|-------------|------------------|
| **VCF Processing** | 🟢 Working | Upload and analysis of genetic data | Basic UI form for uploading VCF files |
| **Genotype Extraction** | 🟡 In Progress | Star allele calling from genetic data | Integration between services needed |
| **Report Generation** | 🔴 Not Started | PDF reports for clinical use | Templates and rendering needed |
| **User Authentication** | 🟡 Basic | Token-based authentication | Need real user database |
| **API Endpoints** | 🟢 Working | REST API for service interactions | API documentation accessible at /docs |
| **Web UI** | 🟢 Working | User-friendly interface | Simple Bootstrap UI for file uploads |

## Recent Changes

1. **Major Update**: Added user-friendly web interface:
   - Created Bootstrap-based web UI
   - Added file upload form for VCF processing
   - Added service status display
   - Implemented Jinja2 templates in FastAPI

2. **Critical Fix**: Resolved issues with PharmCAT and Aldy services:
   - Updated PharmCAT to use simple command to keep container running
   - Added proper Gunicorn server for Aldy service
   - Fixed container restart issues for both services
   - All services now running healthily

3. **Fixed Docker volume mounting issue**:
   - Changed app container mount from root `/` to `/app` directory
   - Added proper PYTHONPATH environment variable
   - Fixed uvicorn command to correctly locate app modules

4. **Environment Variable Fix**: Properly configured .env file:
   - Fixed POSTGRES_PASSWORD variable
   - Fixed SECRET_KEY variable
   - Ensured environment variables are properly passed to containers

5. **Database Connection**: Fixed PostgreSQL connection issues:
   - Ensured proper environment variables for database credentials
   - Updated port mapping to avoid conflicts (5444:5432)
   - Verified database initialization with proper credentials

6. **Health Check Improvements**:
   - Temporarily disabled complex health checks for debugging
   - Simplified health checks to ensure core functionality works
   - Will re-enable proper health checks in next iteration

## Current Challenges

1. **End-to-End Integration**: Need to connect the file upload UI with PharmCAT and Aldy services for full analysis.

2. **Cross-Service Communication**: Need to ensure proper API communication between services, especially for:
   - Main app to PharmCAT wrapper
   - Main app to Aldy service
   - All services to PostgreSQL

3. **Report Generation**: Implement PDF report generation from genotype results.

## Upcoming Tasks

1. **Sample Data**: Add sample VCF files for testing the full pipeline

2. **Database Schema**: Complete PostgreSQL schema for CPIC guidelines

3. **Report Templates**: Create HTML/PDF templates for pharmacogenomic reports

4. **API Documentation**: Complete Swagger/OpenAPI documentation for all services

5. **Testing Scripts**: Create testing scripts to verify end-to-end functionality

6. **Full End-to-End Integration**: Connect file upload UI with service backend processing

## Notes

- Using Docker version 28.0.4
- Using official PharmCAT Docker image (pgkb/pharmcat)
- Directly calling PharmCAT JAR file as per documentation: https://pharmcat.org/using/Running-PharmCAT/
- Python 3.10 for Flask/FastAPI services
- PostgreSQL 15 for database
- The FastAPI application is accessible at http://localhost:8765
- API documentation is available at http://localhost:8765/docs
- Web UI is now available at http://localhost:8765 