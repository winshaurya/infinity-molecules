# Use conda-forge with RDKit pre-installed for faster builds
FROM condaforge/mambaforge:latest

# Set working directory
WORKDIR /app

# Create conda environment with RDKit
RUN mamba create -n rdkit-env python=3.11 rdkit numpy pandas -c conda-forge -y && \
    conda clean --all -y

# Activate environment
ENV PATH="/opt/conda/envs/rdkit-env/bin:$PATH"

# Copy requirements first for better caching
COPY requirements.txt .

# Install remaining Python dependencies (excluding RDKit which is already installed)
RUN pip install --no-cache-dir fastapi uvicorn pydantic supabase python-dotenv networkx psutil

# Copy the rest of the application
COPY . .

# Expose port
EXPOSE 8000

# Run the application
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
