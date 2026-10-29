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
    import re
    from pathlib import Path

    # Recursively find all HTML files staged by Nextflow
    html_files = [f for f in Path('.').rglob('*.html') if f.name != 'index.html']

    # PhylogeneticTreeReport_<subtype>[_<SEGMENT>].html: one whole-genome tree (--phylo-whole-genome) and/or one
    # segment tree per built segment (--phylo-segment-trees), for each --phyloSubtype (e.g. H3N2, H1N1pdm09). The
    # legacy, pre-subtype name PhylogeneticTreeReport.html (older runs' output folders) still matches, with an
    # empty subtype/tree (shown as a single "Phylogenetic Tree" entry, same as before this module built several).
    PHYLO_REPORT_RE = re.compile(r'^phylogenetictreereport(?:_(?P<tag>.+))?\$', re.IGNORECASE)
    # Keep in sync with subworkflows/Phylogenetics.nf's phyloKnownSubtypes() (Groovy/Python, so not shareable
    # directly): both lists gate which subtypes the pipeline and this dashboard recognise.
    PHYLO_KNOWN_SUBTYPES = ["H1N1pdm09", "H3N2"]  # neither is a prefix of the other's "_SEGMENT" suffix

    def split_subtype_tree(tag):
        # "H3N2" (whole genome) -> ("H3N2", "Whole genome"); "H3N2_HA" (segment tree) -> ("H3N2", "HA")
        for s in PHYLO_KNOWN_SUBTYPES:
            if tag == s:
                return s, "Whole genome"
            if tag.startswith(s + "_"):
                return s, tag[len(s) + 1:]
        return tag, ""

    report_catalog = []
    for file_path in html_files:
        clean_name = file_path.stem.replace('_', ' ')

        category = "General Analysis"
        subcategory = ""
        evo_key = ""
        subtype = ""
        tree = ""

        # Interactive tree(s) of the optional phylogenetics module get their own section, one per subtype x tree
        phylo_match = PHYLO_REPORT_RE.match(file_path.stem)
        if phylo_match:
            category = "Phylogenetics"
            tag = phylo_match.group("tag") or ""
            subtype, tree = split_subtype_tree(tag) if tag else ("", "")
            clean_name = f"Phylogenetic Tree {subtype} {tree}".strip()

        # Strictly ensure the filename starts with 'evolution' to exclude clade reports
        if clean_name.lower().startswith('evolution'):
            category = "Frequency Evolution"
            parts = clean_name.split(' ')

            # Extract the protein segment (e.g., HA1, PB1) to create the nested folder
            if len(parts) >= 2:
                subcategory = parts[1].upper()

            # Raw key exactly as emitted by DateGraphicReport.nf's {plot_name}, e.g.
            # "evolution_HA1_A(H3N2).html" -> "HA1_A(H3N2)". This must match the EVO_KEY
            # computed in MutationsGraphicReport.nf / InteractiveMutationsTable.nf verbatim.
            if file_path.stem.lower().startswith('evolution_'):
                evo_key = file_path.stem[len('evolution_'):]

        content = file_path.read_bytes()
        b64_content = base64.b64encode(content).decode('utf-8')

        report_catalog.append({
            "title": clean_name,
            "category": category,
            "subcategory": subcategory,
            "evo_key": evo_key,
            "subtype": subtype,
            "tree": tree,
            "b64": b64_content
        })

    catalog_json = json.dumps(report_catalog)

    segment_mapping = {
        'PB2': 1,
        'PB1': 2,
        'PB1-F2': 2,
        'PA': 3,
        'PA-X': 3,
        'HA1': 4,
        'HA2': 4,
        'NP': 5,
        'NA': 6,
        'M1': 7,
        'M2': 7,
        'NS1': 8,
        'NS2': 8,
    }

    segment_order_json = json.dumps(segment_mapping)

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
            .report-link.loading {{ color: #999; font-style: italic; }}
            
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
            const segmentOrder = {segment_order_json};
            const listContainer = document.getElementById('report-list');
            const iframe = document.getElementById('report-frame');

            // Cache blob URLs 
            function getReportUrl(item) {{
                if (item._blobUrl) return item._blobUrl;

                const byteChars = atob(item.b64);
                const len = byteChars.length;
                const bytes = new Uint8Array(len);
                for (let i = 0; i < len; i++) {{
                    bytes[i] = byteChars.charCodeAt(i);
                }}
                const blob = new Blob([bytes], {{ type: 'text/html' }});
                item._blobUrl = URL.createObjectURL(blob);
                return item._blobUrl;
            }}

            function openReport(item, linkEl) {{
                if (linkEl) {{
                    const originalText = linkEl.textContent;
                    linkEl.classList.add('loading');
                    // Decoding tens of MB of base64 can take a beat; defer to
                    // the next frame so the "loading" state actually paints.
                    requestAnimationFrame(() => {{
                        iframe.src = getReportUrl(item);
                        linkEl.classList.remove('loading');
                    }});
                }} else {{
                    iframe.src = getReportUrl(item);
                }}
            }}

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

                    // Phylogenetics: one report per --phyloSubtype x tree (whole genome and/or one per built
                    // segment). Two dropdowns switch between them instead of listing them as separate links.
                    if (category === 'Phylogenetics' && data.items.length) {{
                        const segTreeOrder = {{ 'PB2': 1, 'PB1': 2, 'PA': 3, 'HA': 4, 'NP': 5, 'NA': 6, 'MP': 7, 'NS': 8 }};
                        const treeRank = t => t === 'Whole genome' ? 0 : (segTreeOrder[t] !== undefined ? segTreeOrder[t] : 999);

                        const selWrapper = document.createElement('div');
                        selWrapper.className = 'report-link';
                        selWrapper.style.cursor = 'default';

                        const subtypeLabel = document.createElement('div');
                        subtypeLabel.textContent = 'Subtype';
                        subtypeLabel.style.cssText = 'font-size:10px;font-weight:bold;color:#666;margin-bottom:4px;';
                        const subtypeSelect = document.createElement('select');
                        subtypeSelect.style.cssText = 'width:100%;padding:6px;border-radius:4px;border:1px solid #ccc;margin-bottom:8px;';

                        const treeLabel = document.createElement('div');
                        treeLabel.textContent = 'Tree';
                        treeLabel.style.cssText = 'font-size:10px;font-weight:bold;color:#666;margin-bottom:4px;';
                        const treeSelect = document.createElement('select');
                        treeSelect.style.cssText = 'width:100%;padding:6px;border-radius:4px;border:1px solid #ccc;';

                        const bySubtype = {{}};
                        data.items.forEach(item => {{
                            const key = item.subtype || item.title;
                            (bySubtype[key] = bySubtype[key] || []).push(item);
                        }});
                        const subtypeOrder = Object.keys(bySubtype).sort();

                        function populateTreeSelect(subtype) {{
                            treeSelect.innerHTML = '';
                            const items = bySubtype[subtype].slice().sort((a, b) => treeRank(a.tree) - treeRank(b.tree));
                            items.forEach((item, idx) => {{
                                const option = document.createElement('option');
                                option.value = String(idx);
                                option.textContent = item.tree || item.title;
                                treeSelect.appendChild(option);
                            }});
                            treeSelect._items = items;
                        }}

                        subtypeOrder.forEach(subtype => {{
                            const option = document.createElement('option');
                            option.value = subtype;
                            option.textContent = subtype;
                            subtypeSelect.appendChild(option);
                        }});
                        populateTreeSelect(subtypeOrder[0]);

                        subtypeSelect.onchange = () => {{
                            populateTreeSelect(subtypeSelect.value);
                            openReport(treeSelect._items[0]);
                        }};
                        treeSelect.onchange = () => openReport(treeSelect._items[parseInt(treeSelect.value, 10)]);

                        selWrapper.appendChild(subtypeLabel);
                        selWrapper.appendChild(subtypeSelect);
                        selWrapper.appendChild(treeLabel);
                        selWrapper.appendChild(treeSelect);
                        catContent.appendChild(selWrapper);
                    }} else {{
                        // Append standalone items directly under the main category
                        data.items.forEach(item => {{
                            const link = document.createElement('a');
                            link.className = 'report-link';
                            link.textContent = item.title;
                            link.onclick = () => openReport(item, link);
                            catContent.appendChild(link);
                        }});
                    }}

                    // Sort subcategories by segment order (genome position); unmapped segments fall to the end, alphabetically among themselves
                    const sortedSubcategories = Object.keys(data.subcategories).sort((a, b) => {{
                        const orderA = segmentOrder[a] !== undefined ? segmentOrder[a] : 999;
                        const orderB = segmentOrder[b] !== undefined ? segmentOrder[b] : 999;
                        if (orderA !== orderB) return orderA - orderB;
                        return a.localeCompare(b);
                    }});
                    
                    for (const subcat of sortedSubcategories) {{
                        const subitems = data.subcategories[subcat].slice().sort((a, b) => {{
                            const numA = parseInt((a.title.match(/H(\\d+)\\s*\$/i) || [])[1] || '999', 10);
                            const numB = parseInt((b.title.match(/H(\\d+)\\s*\$/i) || [])[1] || '999', 10);
                            return numA - numB;
                        }});
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
                            link.onclick = () => openReport(item, link);
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

            // Cross-report navigation broker. MutationsReport.html and MutationsTable.html
            // run inside blob: URI iframes with an opaque-ish origin, so they cannot navigate
            // each other directly - they postMessage up to us instead, and we mediate.
            function openEvolutionReport(evoKey, mutation) {{
                var match = catalog.find(function(item) {{ return item.evo_key === evoKey; }});
                if (!match) {{
                    console.warn('No Frequency Evolution report found for key:', evoKey);
                    return;
                }}

                // Forward the target mutation into the new report only after it has
                // actually finished loading (its own Plotly render happens on load too).
                iframe.onload = function() {{
                    iframe.contentWindow.postMessage(
                        {{ type: 'selectMutation', mutation: mutation }}, '*'
                    );
                    iframe.onload = null;
                }};
                iframe.src = getReportUrl(match);
            }}

            window.addEventListener('message', function(event) {{
                if (!event.data || event.data.type !== 'openEvolution') return;
                openEvolutionReport(event.data.evoKey, event.data.mutation);
            }});

            renderList();
            if (catalog.length > 0) iframe.src = getReportUrl(catalog[0]);
        </script>
    </body>
    </html>'''

    with open('index.html', 'w', encoding='utf-8') as f:
        f.write(html_content)
    """
}
