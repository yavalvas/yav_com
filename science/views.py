from django.contrib.auth.models import User, Group

from django.http import HttpResponse
from django.template import RequestContext, loader
#from django.views.decorators.clickjacking import xframe_options_exempt
#from djangosecure.decorators import frame_deny_exempt


#@xframe_options_exempt
def index(request):
    template = loader.get_template('library/index.html')
    context = RequestContext(request, {})
    return HttpResponse(template.render(context))

#@xframe_options_exempt
def ya_maps(request):
    template = loader.get_template('library/ya_map.html')
    context = RequestContext(request, {})
    return HttpResponse(template.render(context))
