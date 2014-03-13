from django.conf.urls import patterns, include, url
from bio.view import *

#from django.contrib import admin
#admin.autodiscover()

urlpatterns = patterns('',
    ("^hello/$",hello),
    ("^date/$",date),
    # Examples:
    # url(r'^$', 'bio.views.home', name='home'),
    # url(r'^blog/', include('blog.urls')),
#    url(r'^admin/', include(admin.site.urls)),
)
