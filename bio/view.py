from django.http import HttpResponse
import datetime

def hello(request):
    return HttpResponse("Hello world")

def date(request):
    now=datetime.datetime.now()
    return HttpResponse("Now is %s" %now)
