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
| **Main App (FastAPI)** | 🟡 Building | Main application handling UI, reports, user auth | Health check issues being resolved |
| **PostgreSQL Database** | 🟢 Working | Stores CPIC guidelines, user data, reports | Need to add seed data scripts |
| **PharmCAT Service** | 🟢 Working | Using official pgkb/pharmcat Docker image | Using CLI mode instead of API |
| **PharmCAT Wrapper** | 🟡 Building | Python wrapper providing REST API for PharmCAT | Now executes PharmCAT container via Docker CLI |
| **Aldy Service** | 🟡 Building | CYP2D6 specific genotype caller | Need to verify functionality |

### 2. Docker Infrastructure

| Component | Status | Description | Issues/Next Steps |
|-----------|--------|-------------|------------------|
| **Docker Compose** | 🟡 In Progress | Multi-container orchestration | Health check configuration optimizations |
| **Network Config** | 🟢 Working | PGX-Network for service communication | - |
| **Volumes** | 🟢 Working | Data persistence for PostgreSQL | Need additional volumes for shared files |
| **Health Checks** | 🟠 Problematic | Container health monitoring | Fixing issues with Flask/Gunicorn containers |

### 3. Features & Capabilities

| Feature | Status | Description | Issues/Next Steps |
|---------|--------|-------------|------------------|
| **VCF Processing** | 🟡 In Progress | Upload and analysis of genetic data | Need to test with sample data |
| **Genotype Extraction** | 🟡 In Progress | Star allele calling from genetic data | Integration between services needed |
| **Report Generation** | 🔴 Not Started | PDF reports for clinical use | Templates and rendering needed |
| **User Authentication** | 🟡 Basic | Token-based authentication | Need real user database |
| **API Endpoints** | 🟡 In Progress | REST API for service interactions | Documentation needed |

## Recent Changes

1. **Critical Fix**: Updated PharmCAT integration:
   - Removed incorrect `serve --port 8080` command from the PharmCAT container
   - Switched from REST API approach to direct CLI execution
   - Added Docker CLI to PharmCAT wrapper to execute commands in the PharmCAT container
   - Added Docker socket mount to the wrapper container

2. **Major Update**: Switched to the official pgkb/pharmcat Docker image:
   - Removed custom PharmCAT Dockerfile in favor of the maintained image
   - Updated PharmCAT wrapper to communicate with PharmCAT
   - Updated volume mappings to match the official image expectations

3. Fixed Java dependency issues in Dockerfiles:
   - Changed from `openjdk-11-jre-headless` to `default-jre`

4. Created PharmCAT service infrastructure:
   - Created Python wrapper with REST API
   - Added health check scripts

5. Fixed container health check issues:
   - Created custom health check scripts
   - Added proper health check endpoints to services
   - Increased timeout periods for slow-starting services

6. Fixed dependency issues:
   - Added proper requirements.txt for Python services
   - Fixed package dependency resolution (netcat-openbsd)

## Current Challenges

1. **Container Health Checks**: Container health checks are failing, particularly for PharmCAT wrapper service. Investigation ongoing.

2. **Cross-Service Communication**: Need to ensure proper API communication between services, especially for:
   - Main app to PharmCAT wrapper
   - Main app to Aldy service
   - All services to PostgreSQL

3. **Environment Configuration**: Need proper .env file configuration for production deployment.

## Upcoming Tasks

1. **Sample Data**: Add sample VCF files for testing the full pipeline

2. **Database Schema**: Complete PostgreSQL schema for CPIC guidelines

3. **Report Templates**: Create HTML/PDF templates for pharmacogenomic reports

4. **API Documentation**: Generate Swagger/OpenAPI documentation for all services

5. **Testing Scripts**: Create testing scripts to verify end-to-end functionality

## Notes

- Using Docker version 28.0.4
- Using official PharmCAT Docker image (pgkb/pharmcat) in CLI mode
- Python 3.10 for Flask/FastAPI services
- PostgreSQL 15 for database 