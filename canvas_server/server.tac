# -*- coding: utf-8 -*-

"""Приложение Twisted."""

# Поприветствуем пользователя

from twisted.internet import reactor



# Инициализируем сервер Twisted

import canvas_server
canvas_server.setUp()

# Приложение Twisted

application = canvas_server.application
