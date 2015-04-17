import os
import sys
sys.path = ['/root/yavalvas_com'] + sys.path
sys.path.append('/root/yavalvas_com')
sys.path.append('/root/yavalvas_com/graphics')
os.environ['DJANGO_SETTINGS_MODULE'] = 'yavalvas.settings'
import django.core.handlers.wsgi
application = django.core.handlers.wsgi.WSGIHandler()
