from django.conf.urls import patterns, include, url
import graphics.views

urlpatterns = patterns('',
    # Examples:
    # url(r'^$', 'matplot_django.views.home', name='home'),
    # url(r'^blog/', include('blog.urls')),
    (r'^simple$', 'graphics.views.simple'),
    (r'^3D_MODEL$', 'graphics.views.animate'),
    (r'^graphic$', 'graphics.views.graphic'),
    (r'^plot_diff$', 'graphics.views.plot_diff'),
    (r'^all_by_speed$', 'graphics.views.all_by_speed')
)
