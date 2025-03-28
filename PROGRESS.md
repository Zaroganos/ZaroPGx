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
| **Main App (FastAPI)** | 🟢 Working | Main application handling UI, reports, user auth | API is accessible at port 8765 |
| **PostgreSQL Database** | 🟢 Working | Stores CPIC guidelines, user data, reports | Need to add seed data scripts |
| **PharmCAT Service** | 🟡 Building | Using official pgkb/pharmcat Docker image | Container restarting, needs further investigation |
| **PharmCAT Wrapper** | 🟢 Working | Python wrapper providing REST API for PharmCAT | API accessible via port 5001 |
| **Aldy Service** | 🟡 Building | CYP2D6 specific genotype caller | Container restarting, needs further investigation |

### 2. Docker Infrastructure

| Component | Status | Description | Issues/Next Steps |
|-----------|--------|-------------|------------------|
| **Docker Compose** | 🟢 Working | Multi-container orchestration | Core services up and running |
| **Network Config** | 🟢 Working | PGX-Network for service communication | Services communicating successfully |
| **Volumes** | 🟢 Working | Data persistence and sharing | Using volume mounts for JAR sharing |
| **Health Checks** | 🟡 In Progress | Container health monitoring | Temporarily disabled for debugging purposes |

### 3. Features & Capabilities

| Feature | Status | Description | Issues/Next Steps |
|---------|--------|-------------|------------------|
| **VCF Processing** | 🟡 In Progress | Upload and analysis of genetic data | Need to test with sample data |
| **Genotype Extraction** | 🟡 In Progress | Star allele calling from genetic data | Integration between services needed |
| **Report Generation** | 🔴 Not Started | PDF reports for clinical use | Templates and rendering needed |
| **User Authentication** | 🟡 Basic | Token-based authentication | Need real user database |
| **API Endpoints** | 🟢 Working | REST API for service interactions | API documentation accessible at /docs |

## Recent Changes

1. **Critical Fix**: Fixed Docker volume mounting issue:
   - Changed app container mount from root `/` to `/app` directory
   - Added proper PYTHONPATH environment variable
   - Fixed uvicorn command to correctly locate app modules

2. **Environment Variable Fix**: Properly configured .env file:
   - Fixed POSTGRES_PASSWORD variable
   - Fixed SECRET_KEY variable
   - Ensured environment variables are properly passed to containers

3. **Database Connection**: Fixed PostgreSQL connection issues:
   - Ensured proper environment variables for database credentials
   - Updated port mapping to avoid conflicts (5444:5432)
   - Verified database initialization with proper credentials

4. **Health Check Improvements**:
   - Temporarily disabled complex health checks for debugging
   - Simplified health checks to ensure core functionality works
   - Will re-enable proper health checks in next iteration

5. **Major Update**: Switched to the official pgkb/pharmcat Docker image:
   - Removed custom PharmCAT Dockerfile in favor of the maintained image
   - Updated PharmCAT wrapper to communicate with PharmCAT
   - Updated volume mappings to match the official image expectations

6. Fixed Java dependency issues in Dockerfiles:
   - Changed from `openjdk-11-jre-headless` to `default-jre`
   - Added Java JRE to PharmCAT wrapper container

## Current Challenges

1. **PharmCAT & Aldy Services**: These two services are still experiencing restart issues and require further investigation.

2. **Cross-Service Communication**: Need to ensure proper API communication between services, especially for:
   - Main app to PharmCAT wrapper
   - Main app to Aldy service
   - All services to PostgreSQL

3. **Health Checks**: Re-implement proper health checks after core functionality is verified.

## Upcoming Tasks

1. **Sample Data**: Add sample VCF files for testing the full pipeline

2. **Database Schema**: Complete PostgreSQL schema for CPIC guidelines

3. **Report Templates**: Create HTML/PDF templates for pharmacogenomic reports

4. **API Documentation**: Complete Swagger/OpenAPI documentation for all services

5. **Testing Scripts**: Create testing scripts to verify end-to-end functionality

6. **Fix Restarting Services**: Investigate and resolve issues with PharmCAT and Aldy containers

## Notes

- Using Docker version 28.0.4
- Using official PharmCAT Docker image (pgkb/pharmcat)
- Directly calling PharmCAT JAR file as per documentation: https://pharmcat.org/using/Running-PharmCAT/
- Python 3.10 for Flask/FastAPI services
- PostgreSQL 15 for database
- The FastAPI application is accessible at http://localhost:8765
- API documentation is available at http://localhost:8765/docs 