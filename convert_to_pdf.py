#!/usr/bin/env python3
"""
Convert MOECAM Week 2 Progress Report from Markdown to PDF
"""

import markdown
from weasyprint import HTML, CSS
from weasyprint.text.fonts import FontConfiguration
import os

def markdown_to_pdf(md_file, pdf_file):
    """Convert markdown file to PDF with professional styling"""

    # Read markdown content
    with open(md_file, 'r', encoding='utf-8') as f:
        md_content = f.read()

    # Convert markdown to HTML
    html_content = markdown.markdown(md_content, extensions=['tables', 'codehilite'])

    # Create complete HTML document with styling
    full_html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <title>MOECAM Week 2 Progress Report</title>
    </head>
    <body>
        {html_content}
    </body>
    </html>
    """

    # Define CSS styling
    css_style = CSS(string="""
        @page {
            size: A4;
            margin: 2cm;
            @top-center {
                content: "MOECAM Project - Week 2 Progress Report";
                font-family: Arial, sans-serif;
                font-size: 10pt;
                color: #666;
            }
            @bottom-center {
                content: "Page " counter(page) " of " counter(pages);
                font-family: Arial, sans-serif;
                font-size: 10pt;
                color: #666;
            }
        }

        body {
            font-family: Arial, sans-serif;
            font-size: 11pt;
            line-height: 1.6;
            color: #333;
            max-width: 100%;
        }

        h1 {
            font-size: 20pt;
            font-weight: bold;
            color: #2c3e50;
            margin-top: 30pt;
            margin-bottom: 20pt;
            border-bottom: 2pt solid #3498db;
            padding-bottom: 10pt;
        }

        h2 {
            font-size: 16pt;
            font-weight: bold;
            color: #34495e;
            margin-top: 25pt;
            margin-bottom: 15pt;
            border-bottom: 1pt solid #bdc3c7;
            padding-bottom: 5pt;
        }

        h3 {
            font-size: 14pt;
            font-weight: bold;
            color: #2c3e50;
            margin-top: 20pt;
            margin-bottom: 10pt;
        }

        h4 {
            font-size: 12pt;
            font-weight: bold;
            color: #34495e;
            margin-top: 15pt;
            margin-bottom: 8pt;
        }

        p {
            margin-bottom: 10pt;
            text-align: justify;
        }

        ul, ol {
            margin-bottom: 10pt;
            padding-left: 20pt;
        }

        li {
            margin-bottom: 5pt;
        }

        code {
            font-family: 'Courier New', monospace;
            font-size: 10pt;
            background-color: #f8f9fa;
            padding: 2pt 4pt;
            border-radius: 3pt;
            border: 1pt solid #e9ecef;
        }

        pre {
            font-family: 'Courier New', monospace;
            font-size: 9pt;
            background-color: #f8f9fa;
            padding: 10pt;
            border-radius: 5pt;
            border: 1pt solid #e9ecef;
            margin: 10pt 0;
            overflow: hidden;
            line-height: 1.4;
        }

        table {
            width: 100%;
            border-collapse: collapse;
            margin: 15pt 0;
            font-size: 10pt;
        }

        th, td {
            border: 1pt solid #bdc3c7;
            padding: 8pt;
            text-align: left;
        }

        th {
            background-color: #ecf0f1;
            font-weight: bold;
            color: #2c3e50;
        }

        tr:nth-child(even) {
            background-color: #f8f9fa;
        }

        strong {
            font-weight: bold;
            color: #2c3e50;
        }

        hr {
            border: none;
            border-top: 2pt solid #bdc3c7;
            margin: 20pt 0;
        }

        blockquote {
            margin: 15pt 0;
            padding: 10pt 20pt;
            background-color: #f8f9fa;
            border-left: 4pt solid #3498db;
            font-style: italic;
        }

        .page-break {
            page-break-before: always;
        }
    """)

    # Create PDF
    font_config = FontConfiguration()
    html_doc = HTML(string=full_html)
    html_doc.write_pdf(pdf_file, stylesheets=[css_style], font_config=font_config)

    print(f"PDF successfully generated: {pdf_file}")

def main():
    """Main function to convert the progress report"""

    # File paths
    current_dir = os.getcwd()
    md_file = os.path.join(current_dir, "MOECAM_Week2_Progress_Report.md")
    pdf_file = os.path.join(current_dir, "MOECAM_Week2_Progress_Report.pdf")

    # Check if markdown file exists
    if not os.path.exists(md_file):
        print(f"Error: Markdown file not found: {md_file}")
        return

    try:
        # Convert markdown to PDF
        markdown_to_pdf(md_file, pdf_file)

        # Verify PDF was created
        if os.path.exists(pdf_file):
            file_size = os.path.getsize(pdf_file)
            print(f"Success! PDF created: {pdf_file}")
            print(f"File size: {file_size:,} bytes")
        else:
            print("Error: PDF file was not created successfully")

    except Exception as e:
        print(f"Error during PDF conversion: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
