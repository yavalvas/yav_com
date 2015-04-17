from django.conf.urls import patterns, url, include
#from graphics import views
from science import views

from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('',
                       url(r'^graphics/', include('graphics.urls')),
                       url(r'^science/', include('science.urls')),
                       url(r'^canvas_game/', include('canvasgame.urls')),
                       url(r'^admin/', include(admin.site.urls)),
                       # url(r'^', include(router.urls)),
                       # url(r'^api-auth/', include(
                       #     'rest_framework.urls', namespace='rest_framework'))
                       )
