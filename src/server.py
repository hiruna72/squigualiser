import html
from http import HTTPStatus
import http.server
import io
import posixpath
import socketserver
import os
import sys
from urllib.parse import urlparse, parse_qs, quote, unquote, urlsplit, urlunsplit
from plot import argparser, run

class SquigHttpRequestHandler(http.server.SimpleHTTPRequestHandler):
    def do_GET(self):
        if self.path == '/home':
            self.path = 'src/server/home.html'

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
    
    def list_directory(self, path):
        """Helper to produce a directory listing (absent index.html).

        Return value is either a file object, or None (indicating an
        error).  In either case, the headers are sent, making the
        interface the same as for send_head().

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
        try:
            displaypath = unquote(self.path,
                                               errors='surrogatepass')
        except UnicodeDecodeError:
            displaypath = unquote(path)
        displaypath = html.escape(displaypath, quote=False)
        enc = sys.getfilesystemencoding()
        title = 'Directory listing for %s' % displaypath
        r.append('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" '
                 '"http://www.w3.org/TR/html4/strict.dtd">')
        r.append('<html>\n<head>')
        r.append('<meta http-equiv="Content-Type" '
                 'content="text/html; charset=%s">' % enc)
        r.append('<link rel="stylesheet" href="/src/server/style.css">')
        r.append('<title>%s</title>\n</head>' % title)
        r.append('<body>\n<b>%s</b>' % title)
        r.append('<hr>\n<ul>')

        parent_path = '/'.join(self.path.split('/')[0:-2])
        if parent_path == '':
            parent_path = '/'
        if self.path != '/':
            r.append('<li><a href="%s">..</a></li>'
                        % (quote(parent_path,
                                            errors='surrogatepass')))
            
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
                r.append('<li><a target="_blank" href="%s">%s</a></li>'
                        % (quote(linkname,
                                            errors='surrogatepass'),
                        html.escape(displayname, quote=False)))
            else:
                r.append('<li><a href="%s">%s</a></li>'
                    % (quote(linkname,
                                        errors='surrogatepass'),
                    html.escape(displayname, quote=False)))
        r.append('</ul>\n<hr>\n</body>\n</html>\n')
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
print(f'Link to access Squigualiser on browser: http://localhost:{PORT}/home')

my_server.serve_forever()