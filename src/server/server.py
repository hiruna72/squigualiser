import html
from http import HTTPStatus
import http.server
import io
import posixpath
import socketserver
import os
from urllib.parse import urlparse, parse_qs, quote, unquote, urlsplit, urlunsplit

import sys
sys.path.append('../')
from plot import argparser, run

class SquigHttpRequestHandler(http.server.SimpleHTTPRequestHandler):
    def do_GET(self):
        parsed_url = urlparse(self.path)

        # for normal GET requests to be served properly after the translate_path override
        if self.path.startswith('/'):
            self.path = "." + self.path

        if parsed_url.path == '/':
            self.path = 'index.html'

        if parsed_url.path == '/dir_listing':
            query_components = parse_qs(parsed_url.query)
            if 'dir_path' in query_components:
                dir_path = query_components["dir_path"][0]
            else:
                dir_path = "./"

            f = self.list_directory(dir_path)
            if f:
                try:
                    self.copyfile(f, self.wfile)
                finally:
                    f.close()

            return

        if parsed_url.path == "/show_file":
            query_components = parse_qs(parsed_url.query)
            if 'file_path' in query_components:
                file_path = query_components["file_path"][0]
                self.path = file_path
                print(self.path)

        return http.server.SimpleHTTPRequestHandler.do_GET(self)

    def do_POST(self):
        if self.path == "/generate_plots":
            content_len = int(self.headers.get('content-length', 0))
            post_data_bytes = self.rfile.read(content_len)
            plot_params = post_data_bytes.decode("UTF-8").split()
            print(plot_params)

            parser = argparser()
            args = parser.parse_args(plot_params)
            run(args) # TODO: get a return value send it in the response

            self.send_response(200)
            self.send_header('Content-type','text/html')
            self.end_headers()
            self.wfile.write(bytes("Success", "utf8"))

        return

    def translate_path(self, path):
        """Overridden method to allow accessing any path in the system.
           Note: this causes GET requests to be served with 

        """
        # abandon query parameters
        path = path.split('?',1)[0]
        path = path.split('#',1)[0]
        # Don't forget explicit trailing slash when normalizing. Issue17324
        trailing_slash = path.rstrip().endswith('/')
        try:
            path = unquote(path, errors='surrogatepass')
        except UnicodeDecodeError:
            path = unquote(path)
        path = posixpath.normpath(path)
        if trailing_slash:
            path += '/'
        return path
    
    def list_directory(self, path):
        """Overridden method to have a custom formatting for the directory listing html page.

        """
        try:
            list = os.listdir(path)
        except OSError:
            self.send_error(
                HTTPStatus.NOT_FOUND,
                "No permission to list directory")
            return None
        list.sort(key=lambda a: a.lower())
        r = []
        displaypath = unquote(path)
        displaypath = html.escape(displaypath, quote=False)
        enc = sys.getfilesystemencoding()
        title = 'Directory: %s' % displaypath
        r.append('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" '
                 '"http://www.w3.org/TR/html4/strict.dtd">')
        r.append('<html>\n<head>')
        r.append('<meta http-equiv="Content-Type" '
                 'content="text/html; charset=%s">' % enc)
        r.append('<link rel="stylesheet" href="style.css">')
        r.append('<title>%s</title>\n</head>' % title)
        r.append('<body>\n<b>%s</b>' % title)
        r.append('<ul>')
        for name in list:
            fullname = os.path.join(path, name)
            
            displayname = linkname = name
            # Append / for directories or @ for symbolic links
            if os.path.isdir(fullname):
                displayname = name + "/"
                linkname = name + "/"
            if os.path.islink(fullname):
                displayname = name + "@"
                # Note: a link to a directory displays with @ and links with /

            if os.path.isfile(fullname) :
                r.append('<li><a target="_blank" href="show_file?file_path=%s">%s</a></li>'
                        % (quote(fullname,
                                            errors='surrogatepass'),
                        html.escape(displayname, quote=False)))
            else:
                r.append('<li><a href="dir_listing?dir_path=%s">%s</a></li>'
                    % (quote(fullname,
                                        errors='surrogatepass'),
                    html.escape(displayname, quote=False)))
                
        r.append('</ul>\n</body>\n</html>\n')
        encoded = '\n'.join(r).encode(enc, 'surrogateescape')
        f = io.BytesIO()
        f.write(encoded)
        f.seek(0)
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-type", "text/html; charset=%s" % enc)
        self.send_header("Content-Length", str(len(encoded)))
        self.end_headers()
        return f

handler_object = SquigHttpRequestHandler

PORT = 8000
socketserver.TCPServer.allow_reuse_address=True
my_server = socketserver.TCPServer(("", PORT), handler_object)
print("Serving at port", PORT)

my_server.serve_forever()