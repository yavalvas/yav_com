# -*- coding: utf-8 -*-
import sys
import traceback
import random
from twisted.web import server, resource
from twisted.internet import reactor
from twisted.internet.defer import maybeDeferred
from twisted.python import log
from twisted.web.server import Site
from twisted.web.static import File
from twisted.application import internet, service
from twisted.web.resource import Resource
from ws_server import BroadcastServerFactory, BroadcastServerProtocol
from autobahn.twisted.websocket import listenWS
web_root = None
application = service.Application("chess")

# from autobahn.websocket import listenWS



#можно ввести параметр debug и будет производиться логгинг


def setUp():
    if len(sys.argv) > 1 and sys.argv[1] == 'debug':
        log.startLogging(sys.stdout)
        debug = True
    else:
        debug = False
    try:
        import autobahn
        import twisted
    except ImportError:
        sys.exit("Install all dependencies")
    root = Resource()
    # root.putChild(constants.WEB_DYNAMIC_BRANCH, resource)
    from autobahn.twisted.resource import WebSocketResource, HTTPChannelHixie76Aware
    from twisted.web.server import Site

    factory = BroadcastServerFactory("ws://127.0.0.1:8888", debug=debug, debugCodePaths=debug)
    #если используется proxy
    #factory.proxy={'host': '192.168.200.105', 'port': '8088'}
    factory.protocol = BroadcastServerProtocol
    factory.setProtocolOptions(allowHixie76=True)
    ws_resource = WebSocketResource(factory)
    root.putChild("ws", ws_resource)
    site = Site(root)
    site.protocol = HTTPChannelHixie76Aware
    listenWS(factory)
    reactor.run()
    #запуск в качестве демона. объект reactor хранится
    # from twisted.scripts.twistd import run
    # run()
if __name__=="__main__":
    setUp()