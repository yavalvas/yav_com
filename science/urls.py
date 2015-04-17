from django.conf.urls import patterns, url, include
from science import views
from django.conf.urls.static import static
from django.conf import settings

# Wire up our API using automatic URL routing.
# Additionally, we include login URLs for the browseable API.
urlpatterns = patterns('',
                       url(r'^$', views.index, name='index'),
                       url(r'^ya_maps$', views.ya_maps, name='ya_maps'),
                       ) + static(settings.MEDIA_URL, document_root = settings.MEDIA_ROOT)
