"""
Django settings for yavalvas project.

For more information on this file, see
https://docs.djangoproject.com/en/1.6/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/1.6/ref/settings/
"""

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
import os
BASE_DIR = os.path.dirname(os.path.dirname(__file__))


# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/1.6/howto/deployment/checklist/

# SECURITY WARNING: keep the secret key used in production secret!

# SECURITY WARNING: don't run with debug turned on in production!
# DEBUG = False
DEBUG = True
#TEMPLATE_DEBUG = True
TEMPLATE_DEBUG = False
TEMPLATE_DIRS = [os.path.join(BASE_DIR, 'templates')]

ALLOWED_HOSTS = ['82.146.53.119']
#ALLOWED_HOSTS =['127.0.0.1']
# Application definition

INSTALLED_APPS = (
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'science',
    # 'rest_framework',
    'graphics',
    'canvasgame',
)

MIDDLEWARE_CLASSES = (
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
)

#X_FRAME_OPTIONS = 'SAMEORIGIN'
X_FRAME_OPTIONS = 'ALLOWALL'
#X_FRAME_OPTIONS = 'DENY'
ROOT_URLCONF = 'yavalvas.urls'

WSGI_APPLICATION = 'yavalvas.wsgi.application'


# Database
# https://docs.djangoproject.com/en/1.6/ref/settings/#databases


DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': 'django_db'
    }
}

# Internationalization
# https://docs.djangoproject.com/en/1.6/topics/i18n/

LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'UTC'

USE_I18N = True

USE_L10N = True

USE_TZ = True


# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/1.6/howto/static-files/

STATIC_URL = '/static/'

STATICFILES_DIRS = (
    os.path.join(BASE_DIR, "static"),
   # os.path.join(BASE_DIR, "science/static")
)
print "BASEDIR", BASE_DIR
#FIXME: change on deploying..or fix with autodetecting path (os.path.join?)
MEDIA_ROOT = BASE_DIR +"/science/media/"
print "MEDIA_ROOT", MEDIA_ROOT
MEDIA_URL = '/media/'
# REST_FRAMEWORK = {
#     'PAGINATE_BY': 25
# }
