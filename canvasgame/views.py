from django.contrib.auth.models import User, Group

from django.http import HttpResponse
from django.template import RequestContext, loader


def canvasgame(request):
    template = loader.get_template('library/game.html')
    context = RequestContext(request, {})
    return HttpResponse(template.render(context))
