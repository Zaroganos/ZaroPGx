# ZaroPGx - Pharmacogenomic Report Generator

A containerized pharmacogenomic report generator that processes genetic data and generates clinical reports based on CPIC guidelines.

## Features

- Process genetic data from various sources (23andMe, VCF)
- Call alleles using PharmCAT and Aldy
- Generate PDF and interactive HTML reports
- HIPAA-compliant data handling
- RESTful API with authentication

## Architecture

The application consists of several containerized services:

1. **PostgreSQL Database**: Stores CPIC guidelines, gene-drug interactions, and patient data
2. **FastAPI Backend**: Handles file uploads, report generation, and API endpoints
3. **PharmCAT Service**: Performs allele calling for 23andMe data
4. **Aldy Service**: Specialized in CYP2D6 calling from VCF files

## Setup

### Prerequisites

- Docker and Docker Compose
- Git

### Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/Zaroganos/ZaroPGx.git
   cd ZaroPGx
   ```

2. Create and configure the `.env` file:
   ```bash
   cp .env.example .env
   # Edit .env file with your desired settings
   ```

3. Build and start the containers:
   ```bash
   docker-compose up -d
   ```

4. Access the API at `http://localhost:8000`
   - API documentation is available at `http://localhost:8000/docs`

## Usage

### Uploading genetic data

```bash
curl -X POST -F 'file=@sample.txt' -F 'patient_identifier=patient123' \
  -H "Authorization: Bearer YOUR_TOKEN" \
  http://localhost:8000/upload/23andme
```

### Generating a report

```bash
curl -X POST -H "Content-Type: application/json" \
  -H "Authorization: Bearer YOUR_TOKEN" \
  -d '{"patient_id":"1","file_id":"1","report_type":"comprehensive"}' \
  http://localhost:8000/reports/generate
```

## Development

### Project Structure

```
ZaroPGx/
├── app/                    # Python application code
│   ├── api/                # API endpoints and database access
│   ├── pharmcat_wrapper/   # PharmCAT integration
│   └── reports/            # Report generation logic
├── db/                     # Database migrations
│   └── migrations/
│       └── cpic/           # CPIC schema and seed data
├── docker/                 # Additional Dockerfiles
└── docker-compose.yml      # Service orchestration
```

### Adding new genes or drugs

1. Edit the SQL seed files in `db/migrations/cpic/`
2. Restart the containers with `docker-compose restart`

## License

[AGPLv3](LICENSE)

## Acknowledgements

- [CPIC](https://cpicpgx.org/) for pharmacogenomic guidelines
- [PharmCAT](https://pharmcat.org/) for allele calling algorithms
- [Aldy](https://github.com/aldy-team/aldy) for CYP2D6 calling

## Dependency Management

### Container Dependencies

The services are organized in a Docker Compose configuration with clear dependency chains:

1. The PostgreSQL database is the foundation service
2. The PharmCAT service runs independently
3. The PharmCAT Wrapper depends on the PharmCAT service
4. The FastAPI application depends on all other services

### Service-Specific Dependencies

Each service manages its dependencies within its own container:

- **PostgreSQL**: Version 15 with initialization scripts in `db/init`
- **PharmCAT**: Java 17 with PharmCAT v3.0.0 JAR file
- **PharmCAT Wrapper**: Python 3.10 with Flask and other requirements in `docker/pharmcat/requirements.txt`
- **Aldy Service**: Python 3.10 with Aldy and necessary ML libraries
- **FastAPI Backend**: Python 3.10 with requirements specified in `requirements.txt`

### Data Sharing

Services share data through:

1. **Shared Volumes**: The `/data` volume is mounted across containers for file exchange
2. **Network Communication**: All services communicate over the `pgx-network` bridge network
3. **Environment Variables**: Container configuration is managed through environment variables

## Getting Started

### Prerequisites

- Docker and Docker Compose
- Git

### Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/Zaroganos/ZaroPGx.git
   cd ZaroPGx
   ```

2. Create a `.env` file with necessary environment variables:
   ```bash
   POSTGRES_PASSWORD=your_secure_password
   SECRET_KEY=your_secret_key
   ```

3. Build and start the containers:
   ```bash
   docker-compose up -d
   ```

4. The services will be available at:
   - FastAPI Application: http://localhost:8000
   - PharmCAT Service: http://localhost:8080
   - PharmCAT Wrapper: http://localhost:5001
   - Aldy Service: http://localhost:5002

## Development

### Adding New Dependencies

1. For Python services, add new requirements to the appropriate requirements.txt file
2. For system dependencies, add them to the appropriate Dockerfile
3. Rebuild the affected containers:
   ```bash
   docker-compose build <service_name>
   docker-compose up -d <service_name>
   ```

### Service Communication

Services communicate with each other using their service names as hostnames:

- Database: `db:5432`
- PharmCAT: `pharmcat:8080`
- PharmCAT Wrapper: `pharmcat-wrapper:5000`
- Aldy: `aldy:5000`

## Contributing

1. Create a feature branch (`git checkout -b feature/amazing-feature`)
2. Commit your changes (`git commit -m 'Add some amazing feature'`)
3. Push to the branch (`git push origin feature/amazing-feature`)
4. Open a Pull Request

## License

This project is licensed under the MIT License - see the LICENSE file for details. 