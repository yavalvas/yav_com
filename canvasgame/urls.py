from django.conf.urls import patterns, url, include
from canvasgame import views
from django.conf.urls.static import static
from django.conf import settings

# Wire up our API using automatic URL routing.
# Additionally, we include login URLs for the browseable API.
urlpatterns = patterns('',
                       url(r'^canvas$', views.canvasgame, name='cangame')
                       )