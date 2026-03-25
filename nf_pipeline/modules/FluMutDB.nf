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
    local_db_path="${params.protocols[params.protocol].resources}/flumut_db.sql"
    
    current_release="v0.0.0"
    if [ -f "\${local_db_path}" ]; then
        cp "\${local_db_path}" ./flumut_db.sql
        raw_ver=\$(grep -oP 'db_version" VALUES\\(\\K[0-9]+,[0-9]+' flumut_db.sql)
        current_release="v\${raw_ver/,/.}"
    fi

    latest_release=\$(curl -fsSL "https://api.github.com/repos/\${public_repo}/releases/latest" | grep -oP '"tag_name": "\\K[^"]+')
    [ -z "\$latest_release" ] && { echo "ERROR: Failed to fetch latest release."; exit 1; }

    cur_clean=\$(echo "\$current_release" | sed 's/^[v.]*//')
    lat_clean=\$(echo "\$latest_release" | sed 's/^[v.]*//')

    newest=\$(printf '%s\\n' "\$cur_clean" "\$lat_clean" | sort -V | tail -1)

    if [[ "\$newest" != "\$cur_clean" ]]; then
        echo "Updating local database (v\${cur_clean}) to latest release (v\${lat_clean})..."
        curl -fsSL "https://raw.githubusercontent.com/\${public_repo}/\${latest_release}/flumut_db.sql" -o "flumut_db.sql" || { echo "ERROR: Download failed."; exit 1; }
        echo "Database updated successfully."
    else
        echo "Local database is up to date (v\${cur_clean})."
    fi

    python3 -c "import sqlite3; conn = sqlite3.connect('flumut_db.sqlite'); conn.executescript(open('flumut_db.sql').read()); conn.close()"
    rm flumut_db.sql
    """
}