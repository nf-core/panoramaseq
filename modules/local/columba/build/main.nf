process COLUMBA_BUILD {
    tag "columba_build"
    label 'process_medium'
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://quay.io/francoaps/columba_vanilla:latest' :
        'quay.io/francoaps/columba_vanilla:latest' }"
    
    input:
    val columba_repo
    
    output:
    path "build_Vanilla", emit: binaries_dir   // entire directory → single channel item, shareable across samples
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    """
    # Check if running in container (pre-built binaries available)
    if [ -f "/opt/columba/build_Vanilla/columba" ] && [ -f "/opt/columba/build_Vanilla/columba_build" ]; then
        echo "Running in Singularity container - using pre-built Columba binaries"
        cp -rL /opt/columba/build_Vanilla ./
        echo "Container binaries copied successfully"
    else
        # Running on SLURM with HPC modules - use external repository
        echo "Running with HPC modules - using external Columba repository"
        if [ -d "${columba_repo}" ]; then
            echo "Using provided Columba repository: ${columba_repo}"
            REPO_DIR="\$(readlink -f ${columba_repo})"
            # Check if already built
            if [ -d "\${REPO_DIR}/build_Vanilla" ] && [ -f "\${REPO_DIR}/build_Vanilla/columba" ] && [ -f "\${REPO_DIR}/build_Vanilla/columba_build" ]; then
                echo "Found existing build_Vanilla directory with binaries - skipping build"
                NEED_BUILD=false
            else
                echo "No existing build found - will build"
                NEED_BUILD=true
            fi
        else
            echo "ERROR: No Columba repository provided and not running in container"
            echo "Please provide columba_repo parameter when using SLURM profile"
            exit 1
        fi
        # Build Columba only if needed
        if [ "\$NEED_BUILD" = "true" ]; then
            echo "Building Columba in \${REPO_DIR}..."
            cd "\${REPO_DIR}"
            bash build_script.sh Vanilla
            cd -
        else
            echo "Using existing Columba binaries from \${REPO_DIR}"
        fi
        # Copy build directory to output location
        echo "Copying binaries from \${REPO_DIR}/build_Vanilla to ./build_Vanilla"
        cp -rL "\${REPO_DIR}/build_Vanilla" ./
    fi
    
    # Verify binaries exist
    ls -lh build_Vanilla/columba*
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        columba: "\$(./build_Vanilla/columba --version 2>&1 | head -1 || echo 'unknown')"
    END_VERSIONS
    """
    
    stub:
    """
    mkdir -p build_Vanilla
    printf '#!/bin/bash\necho "columba stub"' > build_Vanilla/columba
    printf '#!/bin/bash\necho "columba_build stub"' > build_Vanilla/columba_build
    chmod +x build_Vanilla/columba build_Vanilla/columba_build
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        columba: "stub-version"
    END_VERSIONS
    """
}
