from django.http import HttpResponse
import datetime

def hello(request):
    return HttpResponse("Hello world")

def date(request):
    return HttpResponse(str(datetime.datetime.now()))
