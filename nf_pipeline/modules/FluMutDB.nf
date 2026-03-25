process FluMutDB {
    errorStrategy 'ignore'
    debug true
    
    input:
    path(inferred_subtypes)

    output:
    path "flumut_db.sqlite"

    script:
    """
    public_repo="izsvenezie-virology/FluMutDB"
    local_db_path="${params.protocols[params.protocol].resources}/flumut_db.sqlite"
    
    # Fetch the latest release tag from the public GitHub repository right away
    latest_release=\$(curl -fsSL "https://api.github.com/repos/\${public_repo}/releases/latest" | grep -oP '"tag_name": "\\K[^"]+')
    [ -z "\$latest_release" ] && { echo "ERROR: Failed to fetch latest release."; exit 1; }
    lat_clean=\$(echo "\$latest_release" | grep -o '[0-9].*')

    needs_update=false
    display_version=""

    # If no local database exists, skip the checks and trigger a direct download
    if [ ! -f "\${local_db_path}" ]; then
        needs_update=true
        display_version="No version detected"
    else
        # The file exists, so bring it in and check its internal version
        cp "\${local_db_path}" ./flumut_db.sqlite
        
        raw_ver=\$(python3 -c "import sqlite3; \\
            conn = sqlite3.connect('flumut_db.sqlite'); \\
            cur = conn.cursor(); \\
            cur.execute('SELECT major, minor FROM db_version'); \\
            row = cur.fetchone(); \\
            print(str(row[0]) + '.' + str(row[1])) if row else print('0.0')")
            
        cur_clean=\$(echo "\${raw_ver}" | grep -o '[0-9].*')
        display_version="v\${cur_clean}"
        
        # Compare the local version against the latest GitHub version
        newest=\$(printf '%s\\n' "\$cur_clean" "\$lat_clean" | sort -V | tail -1)
        if [[ "\$newest" != "\$cur_clean" ]]; then
            needs_update=true
        fi
    fi

    # Execute the update or compilation only if the flag was set
    if [ "\$needs_update" = true ]; then
        echo "FluMutDB: Updating local database (\$display_version) to latest release (v\${lat_clean})..."
        
        curl -fsSL "https://raw.githubusercontent.com/\${public_repo}/\${latest_release}/flumut_db.sql" -o "flumut_db.sql" || { echo "ERROR: Download failed."; exit 1; }
        
        rm -f flumut_db.sqlite
        python3 -c "import sqlite3; \\
            conn = sqlite3.connect('flumut_db.sqlite'); \\
            conn.executescript(open('flumut_db.sql').read()); \\
            conn.close()"
            
        rm -f flumut_db.sql
        echo "FluMutDB: Database updated successfully."
    else
        echo "FluMutDB: Local database is up to date (\$display_version). No compilation needed."
    fi
    """
}