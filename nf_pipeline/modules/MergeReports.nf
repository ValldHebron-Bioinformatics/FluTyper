process MergeReports {
    input:
    path reports

    output:
    path "index.html", emit: index

    script:
    """
    #!/usr/bin/env python3
    import base64
    import json
    from pathlib import Path

    # Recursively find all HTML files staged by Nextflow
    html_files = [f for f in Path('.').rglob('*.html') if f.name != 'index.html']

    report_catalog = []
    for file_path in html_files:
        clean_name = file_path.stem.replace('_', ' ')
        
        category = "General Analysis"
        subcategory = ""
        
        # Strictly ensure the filename starts with 'evolution' to exclude clade reports
        if clean_name.lower().startswith('evolution'):
            category = "Frequency Evolution"
            parts = clean_name.split(' ')
            
            # Extract the protein segment (e.g., HA1, PB1) to create the nested folder
            if len(parts) >= 2:
                subcategory = parts[1].upper()
            
        content = file_path.read_bytes()
        b64_content = base64.b64encode(content).decode('utf-8')
        
        report_catalog.append({
            "title": clean_name,
            "category": category,
            "subcategory": subcategory,
            "b64": b64_content
        })

    catalog_json = json.dumps(report_catalog)

    html_content = f'''<!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <title>FluTyper Interactive Dashboard</title>
        <style>
            body {{ font-family: system-ui, sans-serif; margin: 0; display: flex; height: 100vh; overflow: hidden; background: #f9f9f9; }}
            #sidebar {{ width: 350px; background: white; border-right: 1px solid #ddd; display: flex; flex-direction: column; }}
            #search-container {{ padding: 20px; border-bottom: 1px solid #eee; }}
            #search-input {{ width: 100%; padding: 10px; border: 1px solid #ccc; border-radius: 4px; box-sizing: border-box; font-size: 16px; }}
            #report-list {{ flex: 1; overflow-y: auto; padding: 0; }}
            
            /* Main Category Header */
            .category-header {{ 
                padding: 12px 20px; font-weight: bold; color: #555; background: #f0f0f0; 
                text-transform: uppercase; font-size: 12px; letter-spacing: 0.5px;
                cursor: pointer; user-select: none; display: flex; align-items: center;
                border-bottom: 1px solid #e0e0e0;
            }}
            .category-header:hover {{ background: #e8e8e8; }}
            .folder-icon {{ display: inline-block; width: 20px; font-size: 10px; transition: transform 0.2s; }}
            
            /* Subcategory Header */
            .subcategory-header {{
                padding: 10px 20px 10px 30px; font-weight: 600; color: #444; background: #f7f7f7;
                font-size: 11px; cursor: pointer; user-select: none; display: flex; align-items: center;
                border-bottom: 1px solid #eee; text-transform: uppercase;
            }}
            .subcategory-header:hover {{ background: #f0f0f0; }}

            /* Child Links */
            .category-items, .subcategory-items {{ display: block; }}
            .report-link {{ 
                display: block; padding: 10px 20px 10px 40px; color: #333; 
                text-decoration: none; border-bottom: 1px solid #f5f5f5; 
                cursor: pointer; transition: background 0.2s; font-size: 14px;
            }}
            .report-link.sub-link {{ padding-left: 50px; }}
            .report-link:hover {{ background: #e9ecef; color: #007BFF; }}
            
            #viewer-container {{ flex: 1; display: flex; flex-direction: column; background: #fff; }}
            iframe {{ flex: 1; width: 100%; border: none; }}
        </style>
    </head>
    <body>
        <div id="sidebar">
            <div id="search-container">
                <input type="text" id="search-input" placeholder="Search reports..." onkeyup="filterReports()">
            </div>
            <div id="report-list"></div>
        </div>
        <div id="viewer-container">
            <iframe id="report-frame" src="about:blank"></iframe>
        </div>

        <script>
            const catalog = {catalog_json};
            const listContainer = document.getElementById('report-list');
            const iframe = document.getElementById('report-frame');

            function renderList(filterText = '') {{
                listContainer.innerHTML = '';
                const term = filterText.toLowerCase();
                
                const grouped = catalog.reduce((acc, item) => {{
                    const titleMatch = item.title.toLowerCase().includes(term);
                    const catMatch = item.category.toLowerCase().includes(term);
                    const subMatch = item.subcategory && item.subcategory.toLowerCase().includes(term);
                    
                    if (titleMatch || catMatch || subMatch) {{
                        if (!acc[item.category]) acc[item.category] = {{ items: [], subcategories: {{}} }};
                        
                        if (item.subcategory) {{
                            if (!acc[item.category].subcategories[item.subcategory]) {{
                                acc[item.category].subcategories[item.subcategory] = [];
                            }}
                            acc[item.category].subcategories[item.subcategory].push(item);
                        }} else {{
                            acc[item.category].items.push(item);
                        }}
                    }}
                    return acc;
                }}, {{}});

                // Explicitly sort categories to force General Analysis to the top
                const sortedCategories = Object.keys(grouped).sort((a, b) => {{
                    if (a === 'General Analysis') return -1;
                    if (b === 'General Analysis') return 1;
                    return a.localeCompare(b);
                }});

                for (const category of sortedCategories) {{
                    const data = grouped[category];
                    const catWrapper = document.createElement('div');
                    
                    const header = document.createElement('div');
                    header.className = 'category-header';
                    header.innerHTML = '<span class="folder-icon">▼</span> ' + category;
                    
                    const catContent = document.createElement('div');
                    catContent.className = 'category-items';
                    
                    header.onclick = () => {{
                        const isHidden = catContent.style.display === 'none';
                        catContent.style.display = isHidden ? 'block' : 'none';
                        header.querySelector('.folder-icon').textContent = isHidden ? '▼' : '▶';
                    }};

                    // Append standalone items directly under the main category
                    data.items.forEach(item => {{
                        const link = document.createElement('a');
                        link.className = 'report-link';
                        link.textContent = item.title;
                        link.onclick = () => {{ iframe.src = "data:text/html;base64," + item.b64; }};
                        catContent.appendChild(link);
                    }});

                    // Sort subcategories alphabetically before appending them
                    const sortedSubcategories = Object.keys(data.subcategories).sort();
                    
                    for (const subcat of sortedSubcategories) {{
                        const subitems = data.subcategories[subcat];
                        const subWrapper = document.createElement('div');
                        
                        const subHeader = document.createElement('div');
                        subHeader.className = 'subcategory-header';
                        subHeader.innerHTML = '<span class="folder-icon">▼</span> ' + subcat;
                        
                        const subContent = document.createElement('div');
                        subContent.className = 'subcategory-items';
                        
                        subHeader.onclick = () => {{
                            const isHidden = subContent.style.display === 'none';
                            subContent.style.display = isHidden ? 'block' : 'none';
                            subHeader.querySelector('.folder-icon').textContent = isHidden ? '▼' : '▶';
                        }};

                        subitems.forEach(item => {{
                            const link = document.createElement('a');
                            link.className = 'report-link sub-link';
                            link.textContent = item.title;
                            link.onclick = () => {{ iframe.src = "data:text/html;base64," + item.b64; }};
                            subContent.appendChild(link);
                        }});

                        subWrapper.appendChild(subHeader);
                        subWrapper.appendChild(subContent);
                        catContent.appendChild(subWrapper);
                    }}

                    catWrapper.appendChild(header);
                    catWrapper.appendChild(catContent);
                    listContainer.appendChild(catWrapper);
                }}
            }}

            function filterReports() {{
                const text = document.getElementById('search-input').value;
                renderList(text);
            }}

            renderList();
            if (catalog.length > 0) iframe.src = "data:text/html;base64," + catalog[0].b64;
        </script>
    </body>
    </html>'''

    with open('index.html', 'w', encoding='utf-8') as f:
        f.write(html_content)
    """
}