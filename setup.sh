#!/bin/bash

# This script performs the one-time setup for the annotation pipeline.
# It MUST be run once AFTER creating and activating the conda environment.

echo "--- Starting one-time setup for the Annotation Pipeline ---"

# 1. Check if we are inside a conda environment
if [ -z "$CONDA_PREFIX" ]; then
    echo "FATAL ERROR: This script must be run from within an activated conda environment."
    echo "Please run 'conda activate assembly_tool_env' first."
    exit 1
fi

echo "Conda environment detected at: $CONDA_PREFIX"

# 2. Find and Configure RepeatMasker
echo -e "\n--- Searching for RepeatMasker library directory... ---"

RM_DIR=$(find "$CONDA_PREFIX" -type d -name "RepeatMasker" 2>/dev/null | head -n 1)

if [ -z "$RM_DIR" ]; then
    echo "FATAL ERROR: Could not find the RepeatMasker library directory."
    exit 1
fi

echo "RepeatMasker found at: $RM_DIR"
echo "--- Configuring RepeatMasker automatically... ---"

cd "$RM_DIR"

if [ -f "RepeatMasker.lib" ]; then
    echo "RepeatMasker appears to be already configured. Skipping."
else
    echo "Running 'perl ./configure'. This may take a few moments..."
    
    # THE FIX IS HERE: Using a "here document" for robust automation.
    perl ./configure <<EOF
2

5
EOF

    # Verify that configuration was successful
    if [ -f "RepeatMasker.lib" ]; then
        echo "--- RepeatMasker configured successfully! ---"
    else
        echo "FATAL ERROR: RepeatMasker configuration failed."
        echo "Please try running 'perl ./configure' manually inside the directory: $RM_DIR"
        exit 1
    fi
fi

echo -e "\n--- Setup complete! You can now run the main assembly_tool.py pipeline. ---"
