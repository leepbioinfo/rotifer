#!/usr/bin/env python3

def display_html_popup_from_file(file_path,open_file =False,  title="Popout Window", use_button=True):
    from IPython.display import display, HTML
    import html
    # Read the content of the HTML or SVG file
    if open_file:
        with open(file_path, 'r', encoding='utf-8') as file:
            html_content = file.read()
    else:
         html_content = file_path

    # Properly escape the content for safe inclusion in the HTML
    escaped_content = html.escape(html_content).replace("\n", "&#10;")

    if use_button:
        # Generate the pop-out script for a new window
        script = f"""
        <script>
            function openPopup() {{
                var newWindow = window.open("", "{title}", "width=800,height=600");
                newWindow.document.write(`
                    <html>
                        <head>
                            <title>{title}</title>
                        </head>
                        <body>
                            {html_content}
                        </body>
                    </html>
                `);
                newWindow.document.close();
            }}
        </script>
        <button onclick="openPopup()">Open sequence alignment  in Popout</button>
        """
    else:
        # Directly display the content embedded in an iframe
        script = f"""
        <iframe style="border:none; width:100%; height:600px;" srcdoc="{escaped_content}"></iframe>
        """

    display(HTML(script))

def is_running_in_jupyter():
    try:
        from IPython import get_ipython
        if 'IPKernelApp' in get_ipython().config:
            return True
        else:
            return False
    except (ImportError, AttributeError):
        return False    
