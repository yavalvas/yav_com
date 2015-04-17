# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('Agg')
from django.shortcuts import render
import random
import django
import datetime
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.dates import DateFormatter
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from django.views.decorators.csrf import csrf_exempt
import matplotlib.pyplot as plt
import pylab
from numpy import pi, sin, cos, tan, arctan, poly1d, sqrt, array, append
from itertools import izip
import matplotlib
import matplotlib.pyplot as plt
import StringIO
import urllib, base64
from django.http import HttpResponse
import json

@csrf_exempt
def simple(request):
    simpleimgdata = StringIO.StringIO()
    simpleimgdata.buf=""
    fig=Figure()
    ax=fig.add_subplot(111)
    x=[]
    y=[]
    now=datetime.datetime.now()
    delta=datetime.timedelta(days=1)
    for i in range(10):
        x.append(now)
        now+=delta
        y.append(random.randint(0, 1000))
    ax.plot_date(x, y, '-')
    ax.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))
    fig.autofmt_xdate()
    plt.savefig(simpleimgdata, format='png')
    simpleimgdata.seek(0)
    uri = 'data:image/png;base64,' + urllib.quote(base64.b64encode(simpleimgdata.buf))
    # element= '<img src = "%s"/>' % uri
    simpleimgdata.truncate()
    simpleimgdata.seek(0)
    simpleimgdata.close()
    response = HttpResponse()
    response.write(uri)
    return response

def animate(request):
        fig = plt.figure()
        ax = p3.Axes3D(fig)
        R1 = 0.085
        R2 = 0.0535
        Co = 299792458.
        Vo=7500
        ax.clear()
        ax.text2D(0.05, 0.95, "Vo=%s m/s"%Vo, transform=ax.transAxes)
        ax.set_xlim(-0.090, 0.090)
        ax.set_ylim(-0.090, 0.090)
        ax.set_zlim(-0.090, 0.090)
        u, v = np.mgrid[0:5 * np.pi:20j, 0:np.pi:30j]
        x = R1 * np.cos(u) * np.sin(v) * np.sqrt(1. - np.power(Vo / Co, 2))
        y = R1 * np.sin(u) * np.sin(v)
        z = R1 * np.cos(v)
        ax.plot_wireframe(x, y, z, color="r")

        x0, y0, z0 = -1, -1, -1
        x1, y1, z1 = 1, 1, 1

        dx = x1 - x0
        dy = y1 - y0
        dz = z1 - z0
        cx, cy, cz = 0, 0, 0
        A = dx * dx + dy * dy + dz * dz
        B = 2 * dx * (x0 - cx) + 2 * dy * (y0 - cy) + 2 * dz * (z0 - cz)
        C = cx * cx + cy * cy + cz * cz + x0 * x0 + y0 * y0 + z0 * z0 - 2 * (cx * x0 + cy * y0 + cz * z0) - R1 * R1
        p0 = np.poly1d([A, B, C])
        t = p0[1]
        u, v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
        x = R2 * np.cos(u) * np.sin(v) * np.sqrt(1. - np.power(Vo / Co, 2))
        y = R2 * np.sin(u) * np.sin(v)
        z = R2 * np.cos(v)
        ax.plot_wireframe(x, y, z, color="blue")
        fig.autofmt_xdate()
        canvas=FigureCanvas(fig)
        response=django.http.HttpResponse(content_type='image/png')
        canvas.print_png(response)
        imgdata = StringIO.StringIO()
        plt.savefig(imgdata, format='png')
        imgdata.seek(0)
        print "Content-type: image/png\n"
        uri = 'data:image/png;base64,' + urllib.quote(base64.b64encode(imgdata.buf))
        # element= '<img src = "%s"/>' % uri
        imgdata.close()
        return HttpResponse(uri)


def solve_disp_eq(betbn, betbt, bet, Znak, c, It, Ia, nb, var):
    """
    Решение дисперсионного уравнения.
    Znak = -1 при преломлении
    Znak = 1 при отражении
    """
    betb = sqrt(betbn ** 2. + betbt ** 2.)
    gamb = 1. / sqrt(1. - betb ** 2.)
    d = c * It / Ia
    Ab = 1. + (nb ** 2. - 1.) * gamb ** 2. * (1. - betbn ** 2.)
    Bb = d ** 2. * (1. - bet ** 2. - (nb ** 2. - 1.) * (gamb * (bet - betbn)) ** 2.)
    Cb = (nb ** 2. - 1.) * gamb ** 2. * d * betbt * (2. - 2. * bet * betbn - d * betbt * (1. - bet ** 2.))
    Qb = Ab - Bb - Cb
    CHb = bet + (nb ** 2. - 1.) * gamb ** 2. * (bet - betbn) * (1. - betbt * d)
    ZNb = 1. - bet ** 2. - (nb ** 2. - 1.) * (gamb * (bet - betbn)) ** 2.
    kbna = Ia * (CHb + Znak * sqrt(Qb)) / (c * ZNb)  # норм.проекция волн.вектора
    kbt = It  # Тангенц.проекция волн.вектора
    iQb = arctan(abs(kbt / kbna))

    wi = kbna * bet * c + Ia  # Частота прел. волны
    ci = wi * cos(iQb) / abs(kbna)  # Скорость света в среде
    if var < 0:
        iQb = -iQb
    # k = kbna / cos(arctan(abs(kbt / kbna)))
    # # ui=betb*c
    # uit = betbt * c
    # uin = betbn * c
    # V = bet * c
    # Ai = -1 / pow(c, 2.) - (pow(nb, 2.) - 1.) * pow(1. - uin / V, 2.) / pow(c, 2.) / (1. - pow(betb, 2.))
    # Bi = -2 * (pow(nb, 2.) - 1.) * (-kbt * uit + Ia * uin / V) * (1. - uin / V) / (pow(c, 2.) * (1. - pow(betb, 2.)))
    # Ci = pow(k, 2.) - (pow(nb, 2.) - 1.) * pow(-kbt * uit + Ia * uin / V, 2.) / pow(c, 2.) / (1. - pow(betb, 2.))
    # '''print "Ai = %s"%Ai
    # print "Bi = %s"%Bi
    # print "Ci = %s"%Ci'''
    # # wi=(-Bi-sqrt(pow(Bi,2)-4*Ai*Ci))/(2*Ai)
    # dispeq = Ai * wi * wi + Bi * wi + Ci
    # '''print "dispeq = %s"%(dispeq,)
    # print "wi= %s"%wi'''
    return (kbna, kbt, iQb, wi, ci)


def chaise_point(x, y, xs1, ts1, ys1, xs2, ts2, ys2):
    """
    Вычисление ближайшего корня к предыдушему
    """
    distance1 = sqrt((xs1 - x) ** 2 + (ys1 - y) ** 2)
    distance2 = sqrt((xs2 - x) ** 2 + (ys2 - y) ** 2)
    if distance1 > distance2:
        # distance = distance2
        result = (xs2, ts2, ys2)
    else:
        # distance = distance1
        result = (xs1, ts1, ys1)
    return result


def solve_disp_eq(betbn, betbt, bet, Znak, c, It, Ia, nb, var):
    """
    Решение дисперсионного уравнения.
    Znak = -1 при преломлении
    Znak = 1 при отражении
    """
    betb = sqrt(betbn ** 2. + betbt ** 2.)
    gamb = 1. / sqrt(1. - betb ** 2.)
    d = c * It / Ia
    Ab = 1. + (nb ** 2. - 1.) * gamb ** 2. * (1. - betbn ** 2.)
    Bb = d ** 2. * (1. - bet ** 2. - (nb ** 2. - 1.) * (gamb * (bet - betbn)) ** 2.)
    Cb = (nb ** 2. - 1.) * gamb ** 2. * d * betbt * (2. - 2. * bet * betbn - d * betbt * (1. - bet ** 2.))
    Qb = Ab - Bb - Cb
    CHb = bet + (nb ** 2. - 1.) * gamb ** 2. * (bet - betbn) * (1. - betbt * d)
    ZNb = 1. - bet ** 2. - (nb ** 2. - 1.) * (gamb * (bet - betbn)) ** 2.
    kbna = Ia * (CHb + Znak * sqrt(Qb)) / (c * ZNb)  # норм.проекция волн.вектора
    kbt = It  # Тангенц.проекция волн.вектора
    iQb = arctan(abs(kbt / kbna))

    wi = kbna * bet * c + Ia  # Частота прел. волны
    ci = wi * cos(iQb) / abs(kbna)  # Скорость света в среде
    if var < 0:
        iQb = -iQb
    # k = kbna / cos(arctan(abs(kbt / kbna)))
    # # ui=betb*c
    # uit = betbt * c
    # uin = betbn * c
    # V = bet * c
    # Ai = -1 / pow(c, 2.) - (pow(nb, 2.) - 1.) * pow(1. - uin / V, 2.) / pow(c, 2.) / (1. - pow(betb, 2.))
    # Bi = -2 * (pow(nb, 2.) - 1.) * (-kbt * uit + Ia * uin / V) * (1. - uin / V) / (pow(c, 2.) * (1. - pow(betb, 2.)))
    # Ci = pow(k, 2.) - (pow(nb, 2.) - 1.) * pow(-kbt * uit + Ia * uin / V, 2.) / pow(c, 2.) / (1. - pow(betb, 2.))
    # '''print "Ai = %s"%Ai
    # print "Bi = %s"%Bi
    # print "Ci = %s"%Ci'''
    # # wi=(-Bi-sqrt(pow(Bi,2)-4*Ai*Ci))/(2*Ai)
    # dispeq = Ai * wi * wi + Bi * wi + Ci
    # '''print "dispeq = %s"%(dispeq,)
    # print "wi= %s"%wi'''
    return (kbna, kbt, iQb, wi, ci)


def chaise_point(x, y, xs1, ts1, ys1, xs2, ts2, ys2):
    """
    Вычисление ближайшего корня к предыдушему
    """
    distance1 = sqrt((xs1 - x) ** 2 + (ys1 - y) ** 2)
    distance2 = sqrt((xs2 - x) ** 2 + (ys2 - y) ** 2)
    if distance1 > distance2:
        # distance = distance2
        result = (xs2, ts2, ys2)
    else:
        # distance = distance1
        result = (xs1, ts1, ys1)
    return result


def calculate_n_lk6(lw):
     A1=2.1391711
     A2=-9.8913489E-3
     A3=8.4704778E-3
     A4=2.8247761E-4
     A5=-1.9072939E-5
     A6=9.3359448E-7
     lw=lw*1000000
     nboo=A1+A2*lw**2+A3/(lw**2)+A4/(lw**4)+A5/(lw**6)+A6/(lw**8);
     nb=sqrt(nboo)
     return nb

def calculate_n_tf5(lw):
     A1=2.9580175
     A2=-8.2686725E-3
     A3=38.383391E-3
     A4=12.219807E-4
     A5=3.1433368E-5
     A6=86.507903E-7
     lw=lw*1000000
     nboo=A1+A2*lw**2+A3/(lw**2)+A4/(lw**4)+A5/(lw**6)+A6/(lw**8);
     nb=sqrt(nboo)
     return nb

def calculate_n1(has_dispersia, w, N1, c = 299792458.):
    if has_dispersia:
        nb = calculate_n_lk6(2*pi*c/w)
    else:
        nb = N1
    return nb

def calculate_n2(has_dispersia, w, N2, c = 299792458.):
    if has_dispersia:
        nb = calculate_n_tf5(2*pi*c/w)
    else:
        nb = N2
    return nb

def make_all(angle_of_start, velocity, R1, R2, N1, N2,  withCoords, withArrows, outside_radius, has_dispersia, full_results="", with_plot_all=False):
    #trg = target
    Vo = velocity
    # while(1):
    #     var = float(raw_input("Input angle [0;67]:"))
    #     if 0 <= var <= 67:
    #         break
    var = 270 - angle_of_start
    theta0 = var * pi / 180  #угол падения на первую грань FIXME: пока работает при углах от 35 до 60
    # R1, R2 = 0.085, 0.0535
    outside_radius = outside_radius
    # x0 = -sqrt(2.)*0.085/2. #TODO: надо будет сделать изменяющимися в зависимости от угла
    # y0 = -sqrt(2.)*0.085
    x1 = -sqrt(2.) * R1 / 2.
    #координаты x1,y2 с учетом, что под 45 в ту точку.TODO: в дальнейшем исправить на варьирующуюся величину
    x1 = cos(pi * angle_of_start / 180) * R1
    y1 = -sqrt(2.) * R1 / 2.
    y1 = sin(pi * angle_of_start / 180) * R1
    # Vo = 7500



    Vd = 0  #проекция скорости вращения линзы на проекцию волнового вектора
    x_center = y_center = 0
    # c = 3.e8 #скорость вадающего в оптическую систему света TODO: по сути она пройдёт атмосфету, будет иной
    c = 299792458.
    #скорость вращения поверхности земли на экваторе
    Vz = 463
    # R1, R2 = 0.085*sqrt(1-Vo**2/c**2), 0.0535*sqrt(1-Vo**2/c**2)

    # theta0 = arcsin((1-(Vo-Vz)**2/c**2)*sin(theta0)/(1+(Vo-Vz)*cos(theta0)/c))
    #учет аберрации
    #пересчет угла в СО, связанную с землей
    # theta0 = arctan(sin(theta0)*sqrt(1-((Vo-Vz)/c)**2)/(cos(theta0)+((Vo-Vz)/c)))

    #trg.write(u"-----------------New calculations for angle %s ---------------------\n" % angle_of_start)

    full_results+=u"-----------------New calculations for angle %s ---------------------\n" % angle_of_start
    #TODO: подобный блок в дальнейшем в функцию запихать
    theta0s = arctan((y_center - y1) / (x_center - x1))
    full_results+=str(theta0s * 180 / pi)+"\n"
    #trg.write(unicode(theta0s * 180 / pi) + "\n")
    alpha0 = theta0s + theta0  #TODO: тут уточнить. Добавить алгоритм, определяющий, когда сложение, когда вычитание
    full_results+="The angle between the direction of propagation of the beam and the axis ox %s"%(alpha0 * 180 / pi, )+"\n"
    #trg.write(u"Угол между направлением распространения луча и осью ox %s\n" % unicode(alpha0 * 180 / pi))
    # A = ((1+(sin(alpha0)-cos(alpha0))/cos(alpha0))**2+1)
    # B = 2*(y1-x1*(1+(sin(alpha0)-cos(alpha0))/cos(alpha0)))*(1+(sin(alpha0)-cos(alpha0))/cos(alpha0))
    # C = -outside_radius**2+(y1-x1*(1+(sin(alpha0)-cos(alpha0))/cos(alpha0)))**2
    # #решаем уравнение
    # p0 = poly1d([A, B, C])
    # print "Корни уравнения \n", p0, "=", 0
    # #trg.write(u"Корни уравнения \n%s = %s\n"%(p0,0))
    # print p0.r
    # #trg.write(unicode(p0.r)+"\n")
    # ts1 = (p0.r[0]-x1)/(c*cos(alpha0))
    # ts2 = (p0.r[1]-x1)/(c*cos(alpha0))
    # ys1 = y1+c*ts1*sin(alpha0)
    # ys2 = y1+c*ts2*sin(alpha0)
    # x0, t0, y0 = chaise_point(x1, y1, p0.r[0], ts1, ys1, p0.r[1], ts2, ys2)
    x0 = x1
    # t0 = (x1-x0)/(c*cos(alpha0))
    y0 = y1 - (outside_radius - R1)
    t0 = (y1 - y0) / c
    # x6 = p5.r[0]
    # t5 = (x6-x5)/(c5*cos(alpha6))
    # y6 = y5+c5*t5*sin(alpha6)
    full_results+="Coordinates of the reflector for constructing of M0: [%s,%s] [m]" % (x0, y0)+"\n"
    #trg.write(u"Координаты точки вне отражателя для построения М0: [%s,%s] [м]\n" % (x0, y0))
    # t0=sqrt((x1-x0)**2+(y1-y0)**2)/c
    full_results+="Beam propagation time from M0 to M1: t0 = %s [sec]" % (t0,)+"\n"
    #trg.write(u"Время распространения луча от M0 к M1: t0 = %s [сек]\n" % (t0,))
    #FIXME: костыль очень грубый
    # if theta0 == 45*pi/180:
    #     x0, y0 = -sqrt(2.)*0.085/2., -sqrt(2.)*0.085



    #Преломление на первой грани, вычисление угла theta1
    # theta1 = snell(1, 1.47, theta0)
    #TODO: радиусы также изменятся относительно земного наблюдателя
    #TODO: также за время распространения луча от земли до спутника сместится сфера-отражатель
    #TODO: кроме того изменится угол падения при пересчете в систему отсчёта, связанную со спутником
    #TODO: изменится частота(длина) волны вследствие эффекта доплера при пересчете в систему отсчета, связанную
    #TODO: со спутником. Ситуация: приёмник отдаляется от источника э-м волны
    lambd = 532e-9
    wo = 2 * pi * c / lambd
    ko = wo / c
    It = ko * sin(theta0)
    # nb = calculate_n1(has_dispersia, c, wo, N1)
    nb = calculate_n1(has_dispersia, wo, N1)

    Znak = -1
    beto = Vo * cos(theta0) / c
    bet = beto
    Ia = wo * (1 - bet * cos(theta0))
    betbn = Vo * cos(theta0) / c
    betbt = Vo * sin(theta0) / c
    values = solve_disp_eq(betbn=betbn, betbt=betbt, bet=bet, Znak=Znak, c=c, It=It, Ia=Ia, nb=nb, var=var)
    kbna, kbt, iQb, w1, c1 = values
    k1 = kbna / cos(arctan(abs(kbt / kbna)))
    theta1 = iQb
    # theta1s = 45*pi/180. #угол между нормалью и осью x
    # alpha1 = theta1 + theta1s
    # print "Угол между направлением распространения луча и осью ox", alpha1*180/pi
    # A = 1.-2.*(Vo+Vd)/c1/cos(alpha1)+((Vo+Vd)**2)/(c1*cos(alpha1))**2+1+2*(sin(alpha1)-cos(alpha1))/cos(alpha1)+\
    #   ((sin(alpha1)-cos(alpha1))/cos(alpha1))**2
    # B = -2*x1*(Vo+Vd)/c1/cos(alpha1)-2*x1*((Vo+Vd)/c1/cos(alpha1))**2+2*(y1-x1-x1*(sin(alpha1)-cos(alpha1))/cos(alpha1))+\
    #   2*(sin(alpha1)-cos(alpha1))/cos(alpha1)*(y1-x1-x1*(sin(alpha1)-cos(alpha1))/cos(alpha1))
    # C = (x1*(Vo+Vd)/(c1*cos(alpha1)))**2+(y1-x1-x1*(sin(alpha1)-cos(alpha1))/cos(alpha1))**2-R2**2
    # #решаем уравнение
    # p1 = poly1d([A, B, C])
    # print "Корни уравнения \n", p1, "=", 0
    # print p1.r
    # ts1 = (p1.r[0]-x1)/(c1*cos(alpha1))
    # ts2 = (p1.r[1]-x1)/(c1*cos(alpha1))
    # ys1 = y1+c1*ts1*sin(alpha1)
    # ys2 = y1+c1*ts2*sin(alpha1)
    # x2, t1, y2 = chaise_point(x1, y1, p1.r[0], ts1, ys1, p1.r[1], ts2, ys2)
    #
    # print "Координаты точки M2: [%s,%s] [м]"%(x2,y2)
    # print "Время распространения луча от M1 к M2: t1 = %s [сек]"%(t1,)
    # p = -(Vo+Vd)*t1 #вспомогат. парамр
    # theta2 = pi/2-abs(arctan((1+(sin(alpha1)-cos(alpha1))/cos(alpha1)+
    #                         (p+x1)/y2)/(1-(1+(sin(alpha1)-cos(alpha1))/cos(alpha1))*((p+x1)/y2))))
    # print "Угол падения theta2 = %s"%(theta2*180./pi,)
    #TODO: добавить решение системы кинематических уравнений движения луча и движения окружности(сферы в будущем??)
    #TODO: незначительно изменится угол падения, смещение по координате x будет приблизительно 1мкр, что тоже даст
    #TODO: определённый вклад
    #TODO: при этом будет по 2 корня обычно. Надо будет находить расстояние между предыдущей точкой и каждой из 2-х найденных
    #TODO: затем определять ближайшую(с кратчайшим расстоянием). Далее находить угол между нормалью в точке и осью ОХ
    #расчет угла падения на вторую грань

    #TODO: подобный блок в дальнейшем в функцию запихать
    theta1s = arctan((y_center - y1) / (x_center - x1))
    full_results+=str(theta1s * 180 / pi)+"\n"
    #trg.write(unicode(theta1s * 180 / pi) + u"\n")
    alpha1 = theta1 + theta1s
    full_results+="The angle between the direction of propagation of the beam and the axis ox %s"%(alpha1 * 180 / pi, )+"\n"
    #trg.write(u"Угол между направлением распространения луча и осью ox %s\n" % unicode(alpha1 * 180 / pi))
    # A = ((1+(sin(alpha1)-cos(alpha1))/cos(alpha1))**2+1)
    # B = 2*(y1-x1*(1+(sin(alpha1)-cos(alpha1))/cos(alpha1)))*(1+(sin(alpha1)-cos(alpha1))/cos(alpha1))
    # C = -R2**2+(y1-x1*(1+(sin(alpha1)-cos(alpha1))/cos(alpha1)))**2
    A = (sin(alpha1) ** 2 / cos(alpha1) ** 2 + (1 - Vo / c1 / cos(alpha1)) ** 2)
    B = (2 * Vo * x1 * (1 - Vo / c1 / cos(alpha1)) / c1 / cos(alpha1) + 2 * (y1 - x1 * sin(alpha1) / cos(alpha1)) * sin(
        alpha1) / cos(alpha1))
    C = (Vo * x1 / c1 / cos(alpha1)) ** 2 - R2 ** 2 + (y1 - x1 * sin(alpha1) / cos(alpha1)) ** 2
    #решаем уравнение
    p1 = poly1d([A, B, C])
    full_results+="".join(["Roots of the equation \n", str(p1), "=", str(0)])+"\n"
    #trg.write(u"Корни уравнения \n%s = %s\n" % (p1, 0))
    full_results+=str(p1.r)+"\n"
    #trg.write(unicode(p1.r) + "\n")
    p_current1 = p1.r[0]
    p_current2 = p1.r[1]
    ts1 = (p_current1 - x1) / (c1 * cos(alpha1))
    ts2 = (p_current2 - x1) / (c1 * cos(alpha1))
    ys1 = y1 + c1 * ts1 * sin(alpha1)
    ys2 = y1 + c1 * ts2 * sin(alpha1)
    x2, t1, y2 = chaise_point(x1, y1, p_current1, ts1, ys1, p_current2, ts2, ys2)

    full_results+="Coordinates of the M2: [%s,%s] [m]" % (x2, y2)+"\n"
    #trg.write(u"Координаты точки M2: [%s,%s] [м]\n" % (x2, y2))
    t1 = sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2) / c1
    full_results+="Beam propagation time from M1 to M2: t1 = %s [sec]" % (t1,)+"\n"
    #trg.write(u"Время распространения луча от M1 к M2: t1 = %s [сек]\n" % (t1,))

    x_center, y_center = Vo * t1, 0
    theta3s = arctan((y_center - y2) / (x_center - x2))
    theta2 = abs(arctan(
        ((y_center - y2) / (Vo * t1 - x2) - tan(alpha1)) / (1 + (tan(alpha1) * (y_center - y2)) / (Vo * t1 - x2))))
    full_results+="The angle between the normal and the axis OX = %s" % (theta3s * 180 / pi,)+"\n"
    #trg.write(u"Угол между нормалью и осью OX = %s\n" % (theta3s * 180 / pi,))


    It = k1 * sin(theta2)
    # nb = N2
    # nb = calculate_n2(has_dispersia, c1, w1, N2)
    nb = calculate_n2(has_dispersia, w1, N2)
    Znak = -1
    beto = Vo * cos(theta2) / c
    bet = beto
    Ia = w1 * (1 - bet * cos(theta2))
    betbn = Vo * cos(theta2) / c
    betbt = Vo * sin(theta2) / c
    values = solve_disp_eq(betbn=betbn, betbt=betbt, bet=bet, Znak=Znak, c=c, It=It, Ia=Ia, nb=nb, var=var)
    kbna, kbt, iQb, w2, c2 = values
    k2 = kbna / cos(arctan(abs(kbt / kbna)))
    theta3 = iQb

    #TODO: подобный блок в дальнейшем в функцию запихать
    # theta3s = arctan((y_center-y2)/(x_center-x2)) #угол между нормалью и ОХ TODO: переместить выше и решать ДУ с использованием
    full_results+=str(theta3s * 180 / pi)+"\n"
    #trg.write(unicode(theta3s * 180 / pi))
    alpha2 = theta3 + theta3s
    full_results+="The angle between the direction of propagation of the beam and the axis ox %s"% (alpha2 * 180 / pi, )+"\n"
    #trg.write(u"Угол между направлением распространения луча и осью ox %s\n" % unicode(alpha2 * 180 / pi))
    # A = ((1+(sin(alpha2)-cos(alpha2))/cos(alpha2))**2+1)
    # B = 2*(y2-x2*(1+(sin(alpha2)-cos(alpha2))/cos(alpha2)))*(1+(sin(alpha2)-cos(alpha2))/cos(alpha2))
    # C = -R2**2+(y2-x2*(1+(sin(alpha2)-cos(alpha2))/cos(alpha2)))**2
    A = (sin(alpha2) ** 2 / cos(alpha2) ** 2 + (1 - Vo / c2 / cos(alpha2)) ** 2)
    B = (2 * (Vo * x2 / c2 / cos(alpha2) - Vo * t1) * (1 - Vo / c2 / cos(alpha2)) + 2 * (
        y2 - x2 * sin(alpha2) / cos(alpha2)) * sin(alpha2) / cos(alpha2))
    C = (Vo * x2 / c2 / cos(alpha2) - Vo * t1) ** 2 - R2 ** 2 + (y2 - x2 * sin(alpha2) / cos(alpha2)) ** 2
    #решаем уравнение
    p2 = poly1d([A, B, C])
    full_results+="".join(["Roots of the equation \n", str(p2), "=", str(0)])+"\n"
    #trg.write(u"Корни уравнения \n%s = %s\n" % (p2, 0))
    full_results+=str(p2.r)+"\n"
    #trg.write(unicode(p2.r) + "\n")
    # ts1 = (p2.r[0]-x2)/(c1*cos(alpha2)) #TODO: определиться здесь с с1 или с2
    # ts2 = (p2.r[1]-x2)/(c1*cos(alpha2))
    # ys1 = y2+c1*ts1*sin(alpha2)
    # ys2 = y2+c1*ts2*sin(alpha2)
    # x3, t2, y3 = chaise_point(x2, y2, p2.r[0], ts1, ys1, p2.r[1], ts2, ys2)
    x3 = p2.r[1]
    t2 = (x3 - x2) / (c2 * cos(alpha2))
    y3 = y2 + c2 * t2 * sin(alpha2)
    full_results+="Coordinates of M3: [%s,%s] [m]" % (x3, y3)+"\n"
    #trg.write(u"Координаты точки M3: [%s,%s] [м]\n" % (x3, y3))
    t2 = sqrt((x3 - x2) ** 2 + (y3 - y2) ** 2) / c2
    full_results+="The propagation time of the beam from M2 to M3: t2 = %s [sec]" % (t2,)+"\n"
    #trg.write(u"Время распространения луча от M2 к M3: t2 = %s [сек]\n" % (t2,))


    #вычисление угла падения на внутреннюю линзу
    # theta4 = theta3
    x_center, y_center = Vo * (t1 + t2), 0
    theta5s = arctan((y_center - y3) / (x_center - x3))
    theta4 = abs(arctan(((y_center - y3) / (Vo * (t1 + t2) - x3) - tan(alpha2)) / (
        1 + (tan(alpha2) * (y_center - y3)) / (Vo * (t1 + t2) - x3))))
    #Преломление на внутренней линзе, вычисление угла theta5
    # theta5 = snell(1.76, 1.47, theta4)
    #Преломление на третьей грани, вычисление угла theta3
    # theta3 = snell(1.47,  1.76, theta2)
    # w1 = wi
    # ko = w1/c1
    #когда использую новую скорость света в среде, возникает ошибка
    It = k2 * sin(theta4)
    # nb = N1
    # nb = calculate_n1(has_dispersia, c2, w2, N1)
    nb = calculate_n1(has_dispersia, w2, N1)
    Znak = -1
    #необходимо тут вычислять угол между нормалью и осью x
    beto = Vo * cos(theta4) / c
    bet = beto
    Ia = w2 * (1 - bet * cos(theta4))
    betbn = Vo * cos(theta4) / c
    betbt = Vo * sin(theta4) / c
    values = solve_disp_eq(betbn=betbn, betbt=betbt, bet=bet, Znak=Znak, c=c, It=It, Ia=Ia, nb=nb, var=var)
    kbna, kbt, iQb, w3, c3 = values
    k3 = kbna / cos(arctan(abs(kbt / kbna)))
    theta5 = iQb

    #TODO: подобный блок в дальнейшем в функцию запихать
    # theta5s = arctan((y_center-y3)/(x_center-x3))
    full_results+=str(theta5s * 180 / pi)+"\n"
    #trg.write(unicode(theta5s * 180 / pi) + "\n")
    alpha3 = theta5s - theta5  #TODO: тут уточнить
    full_results+="The angle between the direction of propagation of the beam and the axis ox %s"%(alpha3 * 180 / pi,)+"\n"
    #trg.write(u"Угол между направлением распространения луча и осью ox %s\n" % unicode(alpha3 * 180 / pi))
    # A = ((1+(sin(alpha3)-cos(alpha3))/cos(alpha3))**2+1)
    # B = 2*(y3-x3*(1+(sin(alpha3)-cos(alpha3))/cos(alpha3)))*(1+(sin(alpha3)-cos(alpha3))/cos(alpha3))
    # C = -R1**2+(y3-x3*(1+(sin(alpha3)-cos(alpha3))/cos(alpha3)))**2
    A = (sin(alpha3) ** 2 / cos(alpha3) ** 2 + (1 - Vo / c3 / cos(alpha3)) ** 2)
    B = (2 * (Vo * x3 / c3 / cos(alpha3) - Vo * (t1 + t2)) * (1 - Vo / c3 / cos(alpha3)) + 2 * (
        y3 - x3 * sin(alpha3) / cos(alpha3)) * sin(alpha3) / cos(alpha3))
    C = (Vo * x3 / c3 / cos(alpha3) - Vo * (t1 + t2)) ** 2 - R1 ** 2 + (y3 - x3 * sin(alpha3) / cos(alpha3)) ** 2
    #решаем уравнение
    p3 = poly1d([A, B, C])
    full_results+="".join(["Roots of the equation \n", str(p3), "=", str(0)])+"\n"
    #trg.write(u"Roots of the equation \n%s = %s\n" % (p3, 0))
    full_results+=str(p3.r)+"\n"
    #trg.write(unicode(p3.r) + "\n")
    ts1 = (p3.r[0] - x3) / (c3 * cos(alpha3))  #TODO: определиться здесь с с1 или с2
    ts2 = (p3.r[1] - x3) / (c3 * cos(alpha3))
    ys1 = y3 + c3 * ts1 * sin(alpha3)
    ys2 = y3 + c3 * ts2 * sin(alpha3)
    x4, t3, y4 = chaise_point(x3, y3, p3.r[0], ts1, ys1, p3.r[1], ts2, ys2)
    # x3 = p2.r[1]
    # t2 = (x3-x2)/(c2*cos(alpha2))
    # y3 = y2+c2*t2*sin(alpha2)
    full_results+="Coordinates of the M4: [%s,%s] [m]" % (x4, y4)+"\n"
    #trg.write(u"Coordinates of the M4: [%s,%s] [m]\n" % (x4, y4))
    t3 = sqrt((x4 - x3) ** 2 + (y4 - y3) ** 2) / c3
    full_results+="Beam propagation time from M3 to M4: t3 = %s [sec]" % (t3,)+"\n"
    #trg.write(u"Время распространения луча от M3 к M4: t3 = %s [сек]\n" % (t3,))

    x_center, y_center = Vo * (t1 + t2 + t3), 0
    theta6s = arctan((y_center - y4) / (x_center - x4))
    theta6 = abs(arctan(((y_center - y4) / (Vo * (t1 + t2 + t3) - x4) - tan(alpha3)) / (
        1 + (tan(alpha3) * (y_center - y4)) / (Vo * (t1 + t2 + t3) - x4))))
    #вычисление угла падения на внешний мениск
    # theta6 = arcsin(R2*sin(pi-theta5)/R1)
    full_results+="Angle of incidence of the external meniscus theta6 = %s"%(theta6 * 180 / pi,)+"\n"
    #trg.write(u"Угол падения на внешний мениск theta6 = %s\n" % unicode(theta6 * 180 / pi))
    #вычисление угла отражения при решении дисп.ур-я
    #когда использую новую скорость света в среде, возникает ошибка
    It = k3 * sin(theta6)
    # nb = N1
    # nb = calculate_n1(has_dispersia, c3, w3, N1)
    nb = calculate_n1(has_dispersia, w3, N1)
    Znak = 1
    #необходимо тут вычислять угол между нормалью и осью x
    beto = Vo * cos(theta6) / c
    bet = beto
    Ia = w3 * (1 - bet * cos(theta6))
    betbn = Vo * cos(theta6) / c
    betbt = Vo * sin(theta6) / c
    values = solve_disp_eq(betbn=betbn, betbt=betbt, bet=bet, Znak=Znak, c=c, It=It, Ia=Ia, nb=nb, var=var)
    kbna, kbt, iQb, w4, c4 = values
    k4 = kbna / cos(arctan(abs(kbt / kbna)))
    theta6 = iQb
    full_results+="The angle of reflection on the external meniscus with the decision disp. equation theta6 = %s"%(theta6 * 180 / pi,)+"\n"
    #trg.write(u"Угол отражения на внешнем мениске с решением дисп. ур-я theta6 =%s\n" % unicode(theta6 * 180 / pi))
    #TODO: подобный блок в дальнейшем в функцию запихать
    # theta6s = arctan((y_center-y4)/(x_center-x4))
    full_results+=str(theta6s * 180 / pi)
    #trg.write(unicode(theta6s * 180 / pi) + "\n")
    alpha4 = theta6s + theta6  #TODO: тут уточнить
    full_results+="The angle between the direction of propagation of the beam and the axis ox %s"%(alpha4 * 180 / pi,)+"\n"
    # #trg.write(u"Угол между направлением распространения луча и осью ox\n"%(alpha4*180/pi,))
    # A = ((1+(sin(alpha4)-cos(alpha4))/cos(alpha4))**2+1)
    # B = 2*(y4-x4*(1+(sin(alpha4)-cos(alpha4))/cos(alpha4)))*(1+(sin(alpha4)-cos(alpha4))/cos(alpha4))
    # C = -R2**2+(y4-x4*(1+(sin(alpha4)-cos(alpha4))/cos(alpha4)))**2
    A = (sin(alpha4) ** 2 / cos(alpha4) ** 2 + (1 - Vo / c4 / cos(alpha4)) ** 2)
    B = (2 * (Vo * x4 / c4 / cos(alpha4) - Vo * (t1 + t2 + t3)) * (1 - Vo / c4 / cos(alpha4)) + 2 * (
        y4 - x4 * sin(alpha4) / cos(alpha4)) * sin(alpha4) / cos(alpha4))
    C = (Vo * x4 / c4 / cos(alpha4) - Vo * (t1 + t2 + t3)) ** 2 - R2 ** 2 + (y4 - x4 * sin(alpha4) / cos(alpha4)) ** 2

    #решаем уравнение
    p4 = poly1d([A, B, C])
    full_results+="".join(["Roots of the equation \n", str(p4), "=", str(0)])+"\n"
    full_results+=str(p4.r)+"\n"
    ts1 = (p4.r[0] - x4) / (c4 * cos(alpha4))
    ts2 = (p4.r[1] - x4) / (c4 * cos(alpha4))
    ys1 = y4 + c4 * ts1 * sin(alpha4)
    ys2 = y4 + c4 * ts2 * sin(alpha4)
    x5, t4, y5 = chaise_point(x4, y4, p4.r[0], ts1, ys1, p4.r[1], ts2, ys2)
    # x3 = p2.r[1]
    # t2 = (x3-x2)/(c2*cos(alpha2))
    # y3 = y2+c2*t2*sin(alpha2)
    full_results+="Coordinates of the M5: [%s,%s] [m]" % (x5, y5)+"\n"
    t4 = sqrt((x5 - x4) ** 2 + (y5 - y4) ** 2) / c4
    full_results+="The propagation time of the beam from the M4 to M5: t4 = %s [sec]" % (t4,)+"\n"



    #вычисление угла падения на внутреннюю линзу
    # alfa = pi - arcsin(R1*sin(theta6)/R2)
    # theta7 = pi - alfa #угол тупой
    x_center, y_center = Vo * (t1 + t2 + t3 + t4), 0
    theta8s = arctan((y_center - y5) / (x_center - x5))
    theta7 = abs(arctan(((y_center - y5) / (Vo * (t1 + t2 + t3 + t4) - x5) - tan(alpha4)) / (
        1 + (tan(alpha4) * (y_center - y5)) / (Vo * (t1 + t2 + t3 + t4) - x5))))


    #вычисление угла преломления на внутренней линзе
    # theta8 = snell(1.47, 1.76, theta7)
    It = k4 * sin(theta7)
    # nb = N2
    # nb = calculate_n2(has_dispersia, c4, w4, N2)
    nb = calculate_n2(has_dispersia, w4, N2)
    Znak = -1
    #необходимо тут вычислять угол между нормалью и осью x
    beto = Vo * cos(theta7) / c
    bet = beto
    Ia = w4 * (1 - bet * cos(theta7))
    betbn = Vo * cos(theta7) / c
    betbt = Vo * sin(theta7) / c
    values = solve_disp_eq(betbn=betbn, betbt=betbt, bet=bet, Znak=Znak, c=c, It=It, Ia=Ia, nb=nb, var=var)
    kbna, kbt, iQb, w5, c5 = values
    k5 = kbna / cos(arctan(abs(kbt / kbna)))
    theta8 = iQb

    #TODO: подобный блок в дальнейшем в функцию запихать
    # theta8s = arctan((y_center-y5)/(x_center-x5))
    full_results+=str(theta8s * 180 / pi)+"\n"
    alpha5 = theta8s + theta8  #TODO: тут уточнить. Добавить алгоритм, определяющий, когда сложение, когда вычитание
    full_results+="The angle between the direction of propagation of the beam and the axis ox %s"%(alpha5 * 180 / pi,)+"\n"
    # A = ((1+(sin(alpha5)-cos(alpha5))/cos(alpha5))**2+1)
    # B = 2*(y5-x5*(1+(sin(alpha5)-cos(alpha5))/cos(alpha5)))*(1+(sin(alpha5)-cos(alpha5))/cos(alpha5))
    # C = -R2**2+(y5-x5*(1+(sin(alpha5)-cos(alpha5))/cos(alpha5)))**2
    A = (sin(alpha5) ** 2 / cos(alpha5) ** 2 + (1 - Vo / c5 / cos(alpha5)) ** 2)
    B = (2 * (Vo * x5 / c5 / cos(alpha5) - Vo * (t1 + t2 + t3 + t4)) * (1 - Vo / c5 / cos(alpha5)) + 2 * (
        y5 - x5 * sin(alpha5) / cos(alpha5)) * sin(alpha5) / cos(alpha5))
    C = (Vo * x5 / c5 / cos(alpha5) - Vo * (t1 + t2 + t3 + t4)) ** 2 - R2 ** 2 + (
        y5 - x5 * sin(alpha5) / cos(alpha5)) ** 2

    #решаем уравнение
    p5 = poly1d([A, B, C])
    full_results+="".join(["Roots of the equation \n", str(p5), "=", str(0)])+"\n"
    full_results+=str(p5.r)+"\n"
    # ts1 = (p5.r[0]-x5)/(c5*cos(alpha5))
    # ts2 = (p5.r[1]-x5)/(c5*cos(alpha5))
    # ys1 = y5+c5*ts1*sin(alpha5)
    # ys2 = y5+c5*ts2*sin(alpha5)
    # x6, t5, y6 = chaise_point(x5, y5, p5.r[0], ts1, ys1, p5.r[1], ts2, ys2)
    # if 33 <= var <= 64: #FIXME: костыль с магическими числами, убрать как-то
    x6 = p5.r[0]
    # else:
    # x6 = p5.r[1]
    t5 = (x6 - x5) / (c5 * cos(alpha5))
    y6 = y5 + c5 * t5 * sin(alpha5)
    full_results+="Coordinates of the M6: [%s,%s] [m]" % (x6, y6)+"\n"
    t5 = sqrt((x6 - x5) ** 2 + (y6 - y5) ** 2) / c5
    full_results+="The propagation time of the beam from the M5 to M6: t5 = %s [sec]" % (t5,)+"\n"






    #вычисление угла падения на внутреннюю линзу
    # theta9 = theta8
    x_center, y_center = Vo * (t1 + t2 + t3 + t4 + t5), 0
    theta10s = arctan((y_center - y6) / (x_center - x6))
    theta9 = abs(arctan(((y_center - y6) / (Vo * (t1 + t2 + t3 + t4 + t5) - x6) - tan(alpha5)) / (
        1 + (tan(alpha5) * (y_center - y6)) / (Vo * (t1 + t2 + t3 + t4 + t5) - x6))))

    #вычисление угла преломления на внутренней линзе
    # theta10 = snell(1.76, 1.47, theta9)
    It = k5 * sin(theta9)
    # nb = N1
    # nb = calculate_n1(has_dispersia, c5, w5, N1)
    nb = calculate_n1(has_dispersia, w5, N1)
    Znak = -1
    #необходимо тут вычислять угол между нормалью и осью x
    beto = Vo * cos(theta9) / c
    bet = beto
    Ia = w5 * (1 - bet * cos(theta9))
    betbn = Vo * cos(theta9) / c
    betbt = Vo * sin(theta9) / c
    values = solve_disp_eq(betbn=betbn, betbt=betbt, bet=bet, Znak=Znak, c=c, It=It, Ia=Ia, nb=nb, var=var)
    kbna, kbt, iQb, w6, c6 = values
    k6 = kbna / cos(arctan(abs(kbt / kbna)))
    theta10 = iQb


    #TODO: подобный блок в дальнейшем в функцию запихать
    # theta10s = arctan((y_center-y6)/(x_center-x6))
    full_results+=str(theta10s * 180 / pi)+"\n"
    alpha6 = theta10s - theta10  #TODO: тут уточнить. Добавить алгоритм, определяющий, когда сложение, когда вычитание
    full_results+="The angle between the direction of propagation of the beam and the axis ox %s"%(alpha6 * 180 / pi,)+"\n"
    # A = ((1+(sin(alpha6)-cos(alpha6))/cos(alpha6))**2+1)
    # B = 2*(y6-x6*(1+(sin(alpha6)-cos(alpha6))/cos(alpha6)))*(1+(sin(alpha6)-cos(alpha6))/cos(alpha6))
    # C = -R1**2+(y6-x6*(1+(sin(alpha6)-cos(alpha6))/cos(alpha6)))**2
    A = (sin(alpha6) ** 2 / cos(alpha6) ** 2 + (1 - Vo / c6 / cos(alpha6)) ** 2)
    B = (2 * (Vo * x6 / c6 / cos(alpha6) - Vo * (t1 + t2 + t3 + t4 + t5)) * (1 - Vo / c6 / cos(alpha6)) + 2 * (
        y6 - x6 * sin(alpha6) / cos(alpha6)) * sin(alpha6) / cos(alpha6))
    C = (Vo * x6 / c6 / cos(alpha6) - Vo * (t1 + t2 + t3 + t4 + t5)) ** 2 - R1 ** 2 + (
        y6 - x6 * sin(alpha6) / cos(alpha6)) ** 2

    #решаем уравнение
    p6 = poly1d([A, B, C])
    full_results+="".join(["Roots of the equation \n", str(p6), "=", str(0)])+"\n"
    full_results+=str(p6.r)+"\n"
    ts1 = (p6.r[0] - x6) / (c6 * cos(alpha6))
    ts2 = (p6.r[1] - x6) / (c6 * cos(alpha6))
    ys1 = y6 + c6 * ts1 * sin(alpha6)
    ys2 = y6 + c6 * ts2 * sin(alpha6)
    x7, t6, y7 = chaise_point(x6, y6, p6.r[0], ts1, ys1, p6.r[1], ts2, ys2)
    # x6 = p5.r[0]
    # t5 = (x6-x5)/(c5*cos(alpha6))
    # y6 = y5+c5*t5*sin(alpha6)
    full_results+="Coords of point M7: [%s,%s] [m]" % (x7, y7)+"\n"
    t6 = sqrt((x7 - x6) ** 2 + (y7 - y6) ** 2) / c6
    full_results+="The propagation time of the beam from the M6 to M7: t6 = %s [sec]" % (t6,)+"\n"



    #вычисление угла падения на внешний мениск
    # theta11 = arcsin(R2*sin(pi-theta10)/R1)
    x_center, y_center = Vo * (t1 + t2 + t3 + t4 + t5 + t6), 0
    theta12s = arctan((y_center - y7) / (x_center - x7))
    theta11 = abs(arctan(((y_center - y7) / (Vo * (t1 + t2 + t3 + t4 + t5 + t6) - x7) - tan(alpha6)) / (
        1 + (tan(alpha6) * (y_center - y7)) / (Vo * (t1 + t2 + t3 + t4 + t5 + t6) - x7))))

    #вычисление угла преломления на внешнем мениске
    # theta12 = snell(1.47, 1., theta11)
    It = k6 * sin(theta11)
    nb = 1
    Znak = -1
    #необходимо тут вычислять угол между нормалью и осью x
    beto = Vo * cos(theta11) / c
    bet = beto
    Ia = w6 * (1 - bet * cos(theta11))
    betbn = Vo * cos(theta11) / c
    betbt = Vo * sin(theta11) / c
    values = solve_disp_eq(betbn=betbn, betbt=betbt, bet=bet, Znak=Znak, c=c, It=It, Ia=Ia, nb=nb, var=var)
    kbna, kbt, iQb, w7, c7 = values
    k7 = kbna / cos(arctan(abs(kbt / kbna)))
    theta12 = iQb
    #учет аберрации
    # theta12 = arctan(sin(theta12)*sqrt(1-((Vo-Vz)/c)**2)/(cos(theta12)+((Vo-Vz)/c)))

    #TODO: подобный блок в дальнейшем в функцию запихать
    # theta12s = arctan((y_center-y7)/(x_center-x7))
    full_results+=str(theta12s * 180 / pi)+"\n"
    alpha7 = theta12s - theta12  #TODO: тут уточнить. Добавить алгоритм, определяющий, когда сложение, когда вычитание
    full_results+="The angle between the direction of propagation of the beam and the axis ox %s"%(alpha7 * 180 / pi,)+"\n"
    A = ((1 + (sin(alpha7) - cos(alpha7)) / cos(alpha7)) ** 2 + 1)
    B = 2 * (y7 - x7 * (1 + (sin(alpha7) - cos(alpha7)) / cos(alpha7))) * (
        1 + (sin(alpha7) - cos(alpha7)) / cos(alpha7))
    C = -outside_radius ** 2 + (y7 - x7 * (1 + (sin(alpha7) - cos(alpha7)) / cos(alpha7))) ** 2
    #решаем уравнение
    p7 = poly1d([A, B, C])
    full_results+="".join(["Roots of the equation \n", str(p7), "=", str(0)])+"\n"
    full_results+=str(p7.r)+"\n"
    ts1 = (p7.r[0] - x7) / (c7 * cos(alpha7))
    ts2 = (p7.r[1] - x7) / (c7 * cos(alpha7))
    ys1 = y7 + c7 * ts1 * sin(alpha7)
    ys2 = y7 + c7 * ts2 * sin(alpha7)
    x8, t7, y8 = chaise_point(x7, y7, p7.r[0], ts1, ys1, p7.r[1], ts2, ys2)
    # x6 = p5.r[0]
    # t5 = (x6-x5)/(c5*cos(alpha6))
    # y6 = y5+c5*t5*sin(alpha6)
    full_results+="Coordinates of the point outside of a reflector for plotting M8: [%s,%s] [m]" % (x8, y8)+"\n"
    t7 = sqrt((x8 - x7) ** 2 + (y8 - y7) ** 2) / c7
    full_results+="The propagation time of the beam from the M7 to M8: t7 = %s [sec]" % (t7,)+"\n"

    for i in xrange(13):
        full_results+="".join(["theta%s" % i, " = ", str(eval("theta%s" % i) * 180 / pi)])+"\n"

    for i in xrange(1, 8):
        full_results+="".join(["c%s" % i, " = ", str(eval("c%s" % i))])+"\n"
        full_results+="".join(["w%s" % i, " = ", str(eval("w%s" % i))])+"\n"
        # print "k%s"%i," = ",eval("k%s"%i)


    # /usr/bin/python /home/pacino/Pycharm_Projects/linze/linze_full.py
    # Координаты точки M1: [-0.0601040764009,-0.0601040764009] [м]
    # угол падения на внешний мениск в градусах theta1 = 45.0
    # угол преломления в градусах theta1 = 28.7523691841
    # Угол между нормалью и осью х = 73.7523691841
    # Корни уравнения
    #        2
    # 12.77 x + 1.003 x + 0.01849 = 0
    # [-0.04891199 -0.02959839]
    # Координаты точки M2: [-0.0489119896281,-0.021699855194] [м]
    # Время распространения луча от M1 к M2: t1 = 1.33339458727e-10 [сек]
    # Угол падения theta2 = 53.9011857132
    # Угол преломления на внутреннем мениске theta3=42.4436963668
    # Угол между нормалью и осью OX = 23.9240543546
    # Угол между направлением распространения луча и OX = 66.3677507213
    # корни уравнения
    #        2
    # 6.223 x + 0.4775 x + 0.008052 = 0
    # [-0.05171194 -0.02502096]
    # Координаты точки M3: [-0.0250209643696,0.0329007503784] [м]
    # Время распространения луча от M2 к M3: t2 = 1.98662394159e-10 [сек]
    # Угол падения theta4 = 46.5167263646
    # Угол преломления на внутренней линзе theta5=60.3097709897
    # Угол между нормалью и осью OX = 52.7441429562
    #
    # Process finished with exit code 0

    # Без дисперсии
    # theta0  =  45.0
    # theta1  =  28.7523691841
    # theta2  =  49.8399026188
    # theta3  =  39.6665533421
    # theta4  =  39.6665533421
    # theta5  =  49.8399026188
    # theta6  =  28.7523691841
    # theta7  =  49.8399026188
    # theta8  =  39.6665533421
    # theta9  =  39.6665533421
    # theta10  =  49.8399026188
    # theta11  =  28.7523691841
    # theta12  =  45.0


    # /usr/bin/python /home/pacino/Pycharm_Projects/linze/calc_linze/linze_main.py
    # theta0  =  60.0
    # theta1  =  36.0955022396
    # theta2  =  69.3912163488
    # theta3  =  51.4236011281
    # theta4  =  51.4236011281
    # theta5  =  69.3912163488
    # theta6  =  36.0955022396
    # theta7  =  69.3912163488
    # theta8  =  51.4236011281
    # theta9  =  51.4236011281
    # theta10  =  69.3912163488
    # theta11  =  36.0955022396
    # theta12  =  60.0
    #
    # Process finished with exit code 0

    # Угол падения на внешний мениск theta6 = 36.3859714639
    # Угол отражения на внешнем мениске с решением дисп. ур-я theta6 = 36.4615277538
    # theta0  =  60.0
    # theta1  =  36.2248888954
    # theta2  =  69.8680207139
    # theta3  =  51.6877075424
    # theta4  =  51.6877075424
    # theta5  =  70.4760572672
    # theta6  =  36.4615277538
    # theta7  =  70.7671683045
    # theta8  =  52.3385126529
    # theta9  =  52.3385126529
    # theta10  =  71.1797563806
    # theta11  =  36.5668638544
    # theta12  =  61.6704582792
    #
    # Process finished with exit code 0

    def plot_all():
        pylab.axes()
        x = [x0, x1, x2, x3, x4, x5, x6, x7, x8]
        y = [y0, y1, y2, y3, y4, y5, y6, y7, y8]
        plt.title('Retroreflective Sphere')
        plt.plot(x, y)
        plt.xlim(-0.1, 0.1)
        plt.ylim(-0.13, 0.1)
        plt.xlabel(u"[m]")
        plt.ylabel(u"[m]")
        point_num = 0
        for i, j in izip(x, y):
            if withCoords:
                plt.annotate("M%s\n" % point_num + "[%2.3e," % i + "%2.3e]" % j, xy=(i, j))
            point_num += 1
        if withArrows:
            for i in xrange(len(x) - 1):
                plt.arrow(x[i],
                          y[i],
                          x[i + 1] - x[i],
                          y[i + 1] - y[i],
                          head_width=0.005,
                          head_length=0.004,
                          fc="k",
                          ec="k",
                          width=0.00003)
    if with_plot_all:
        plot_all()
    x = array([x0, x1, x2, x3, x4, x5, x6, x7, x8])
    y = array([y0, y1, y2, y3, y4, y5, y6, y7, y8])
    list_of_tuples_of_points = zip(x, y)
    return (x, y), theta12 * 180. / pi, theta1s, theta12s, full_results
    # return (x, y), theta12 * 180. / pi, theta1s, theta12s



def draw_dtheta(withSharp=None, Vo=7500, R1=0.085, R2=0.0535, N1=1.47290, N2=1.76470, withCoords=False, withArrows=False, outside_radius=0.12, has_dispersia=False, full_results=""):
    """
    Redraws the figure
    """
    input_angles = array([269.999, 269.5, 269., 268.5, 268., 267.5,
                    267, 266.5, 266, 265.5, 265, 264.5, 264, 263.5, 263,
                    262.5, 262, 260, 259, 258, 257])
    inp_angles = array([270 - i for i in input_angles])
    output_angles_7500 = array([])
    theta1_angles_7500 = array([])
    list_of_theta12s_7500 = array([])
    full_string=full_results
    for data in input_angles:
        angles = make_all(data, Vo, R1, R2, N1, N2, withCoords, withArrows, outside_radius, has_dispersia, "", with_plot_all=False)
        output_angles_7500 = append(output_angles_7500, array([angles[1]]))
        theta1_angles_7500 = append(theta1_angles_7500, array([angles[2]]))
        list_of_theta12s_7500 = append(list_of_theta12s_7500, array([angles[3] * 180. / pi]))
        full_string+=angles[4]
    list_of_output_sum_7500 = array([abs(k) + l for k, l in izip(list_of_theta12s_7500, output_angles_7500)])
    # theta_graduses_7500 = array([i * 180 / pi for i in theta1_angles_7500])
    # main_dtheta_7500 = array([(k + t) for k, t in izip(inp_angles, theta_graduses_7500)])
    output_angles_0 = array([])
    Vo = 0
    list_of_theta12s = array([])
    for data in input_angles:
        angles = make_all(data, Vo, R1, R2, N1, N2, withCoords, withArrows, outside_radius, has_dispersia, "", with_plot_all=False)
        output_angles_0 = append(output_angles_0, array([angles[1]]))
        list_of_theta12s = append(list_of_theta12s, array([angles[3] * 180. / pi]))
        full_string+=angles[4]
    list_of_output_sum_0 = array([abs(k) + l for k, l in izip(list_of_theta12s, output_angles_0)])
    # list_diff_out_0 = array([k - m for k, m in izip(output_angles_0, inp_angles)])
    # list_diff_out_7500 = array([k - m for k, m in izip(output_angles_7500, inp_angles)])
    diff_angles_by_velocity = array([l - m for l, m in izip(output_angles_7500, output_angles_0)])
    # list_of_angles_by_velocity = array([tan(alpha * 180. / pi) * 835000 for alpha in diff_angles_by_velocity])
    list_of_output_dl_0 = array([tan((i - 90) * pi / 180.) * 835000 for i in list_of_output_sum_0])
    list_of_output_dl_7500 = array([tan((i - 90) * pi / 180.) * 835000 for i in list_of_output_sum_7500])
    #Далее построение разницы. Построение при 0 и при 7500 закомментировал
    main_dl_diff = array([i - k for i, k in izip(list_of_output_dl_7500, list_of_output_dl_0)])
    plt.plot(inp_angles, main_dl_diff, linestyle=u'solid', c='y')
    plt.xlabel(u"Input angle [gradus]")
    plt.ylabel(u"dL=dl_7500-dl_0")
    return full_string

def str_to_bool(spec_string):
    return True if spec_string=="true" else False
@csrf_exempt
def graphic(request):
    imgdata = StringIO.StringIO()
    imgdata.truncate(0)
    imgdata.buf=""
    if request.method == 'GET':
        try:
            list_of_angles = map(float, request.GET.get("angles").split(" "))
        except:
            return HttpResponse(json.dumps({"response_error": "wrong_angle_type"}), content_type="application/json")
        for item in list_of_angles:
            if 325<item or item<215:
                return HttpResponse(json.dumps({"response_error": "wrong_angle_value"}), content_type="application/json")
        Vo = float(request.GET.get("Vo"))
        R1 = float(request.GET.get("R1"))
        R2 = float(request.GET.get("R2"))
        N1 = float(request.GET.get("N1"))
        N2 = float(request.GET.get("N2"))
        if R1<R2:
            return HttpResponse(json.dumps({"response_error": "r1_lower_r2"}), content_type="application/json")
        if Vo>3000000 or Vo<-3000000:
            return HttpResponse(json.dumps({"response_error": "wrong_vo"}), content_type="application/json")
        if R1<0 or R2<0 or N1<0 or N2<0:
            return HttpResponse(json.dumps({"response_error": "wrong_value"}), content_type="application/json")
        outside_radius = float(request.GET.get("R_out"))
        has_dispersia = str_to_bool(str(request.GET.get("hasDispersia")))
        withCoords = str_to_bool(str(request.GET.get("withCoords")))
        print "withCoords", type(withCoords), withCoords
        withSharp = str_to_bool(str(request.GET.get("withSharp")))
        if withSharp:
            plt.grid()
        withArrows = str_to_bool(str(request.GET.get("withArrows")))
        full_string=""
        for angle in list_of_angles:
            full_results=make_all(angle, Vo, R1, R2, N1, N2, withCoords, withArrows, outside_radius, has_dispersia, "", with_plot_all=True)[4]
            full_string+=full_results
        circle1 = plt.Circle((0, 0), R1, linewidth=2, facecolor="w")
        circle2 = plt.Circle((0, 0), R2, linewidth=2, facecolor="w")
        fig = plt.gcf()
        if has_dispersia:
            fig.gca().axes.set_title(u"С дисп. $\mathcal{N_1}=\sqrt{{A_1}+{A_2}lw^2+\\frac{A_3}{lw^2}+\\frac{A_4}{lw^4}+\\frac{A_5}{lw^6}+\\frac{A_6}{lw^6}}; \mathcal{N_2}=\sqrt{{A_1}+{A_2}lw^2+\\frac{A_3}{lw^2}+\\frac{A_4}{lw^4}+\\frac{A_5}{lw^6}+\\frac{A_6}{lw^8}}$")
        else:
            fig.gca().axes.set_title(u"Без дисперсии $\mathcal{N_1}=%s, \mathcal{N_2}=%s$"%(N1, N2))
        fig.gca().add_artist(circle1)
        fig.gca().add_artist(circle2)
        fig.autofmt_xdate()
        plt.savefig(imgdata, format='png')
        imgdata.seek(0)
        print "Content-type: image/png\n"
        uri = 'data:image/png;base64,' + urllib.quote(base64.b64encode(imgdata.buf))
        imgdata.truncate()
        imgdata.seek(0)
        imgdata.close()
        # response = HttpResponse()
        # response.write(uri)
        response_data={
            "main_image": uri,
            "full_text": full_string
        }
        return HttpResponse(json.dumps(response_data), content_type="application/json")
        
@csrf_exempt
def plot_diff(request):
    imgdata = StringIO.StringIO()
    imgdata.truncate(0)
    imgdata.buf=""
    if request.method == 'GET':
        withSharp = str(request.GET.get("withSharp"))
        if withSharp=="true":
            plt.grid()
        full_results=""
        full_results = draw_dtheta(withSharp=None, Vo=7500, R1=0.085, R2=0.0535, N1=1.47290, N2=1.76470, withCoords=False, withArrows=False, outside_radius=0.12, has_dispersia=False, full_results=full_results)
        fig = plt.gcf()
        fig.autofmt_xdate()
        plt.savefig(imgdata, format='png')
        imgdata.seek(0)
        print "Content-type: image/png\n"
        uri = 'data:image/png;base64,' + urllib.quote(base64.b64encode(imgdata.buf))
        imgdata.truncate()
        imgdata.seek(0)
        imgdata.close()
        # response = HttpResponse()
        # response.write(uri)
        response_data={
            "main_image": uri,
            "full_text": full_results
        }
        return HttpResponse(json.dumps(response_data), content_type="application/json")



def _calculate_dl_for_angle_and_plot(input_angles, Vo, has_dispersia, full_results):
    output_angles = array([])
    list_of_theta12s = array([])
    R1 = 0.085
    R2 = 0.0535
    N1 = 1.4729
    N2 = 1.7647
    outside_radius = 0.12
    for data in input_angles:
        angles = make_all(data, Vo, R1, R2, N1, N2, False, False, outside_radius, has_dispersia, "", with_plot_all=False)
        output_angles = append(output_angles, array([angles[1]]))
        list_of_theta12s = append(list_of_theta12s, array([angles[3] * 180. / pi]))
        full_results+=angles[4]
    list_of_output_sum = array(
        [abs(k) - l if l < 0 else abs(k) + l for k, l in izip(list_of_theta12s, output_angles)])
    list_of_output_dl = array([tan((i - 90) * pi / 180.) * 835000 for i in list_of_output_sum])
    return list_of_output_dl, full_results


def draw_all_by_speed(input_angles, speed_data, has_dispersia, full_results=""):
    inp_angles = array([270 - i for i in input_angles])
    for Vo in speed_data:
        list_of_output_dl, full_results = _calculate_dl_for_angle_and_plot(input_angles, Vo, has_dispersia, full_results)
        plt.plot(inp_angles, list_of_output_dl)
    plt.xlabel(u"Input angle [gradus]")
    plt.ylabel(u"dl")
    return full_results

@csrf_exempt
def all_by_speed(request):
    imgdata = StringIO.StringIO()
    imgdata.truncate(0)
    imgdata.buf=""
    if request.method == 'GET':
        withSharp = str(request.GET.get("withSharp"))
        if withSharp=="true":
            plt.grid()
        input_angles = map(float, request.GET.get("angles").split(" "))
        speed_data = map(float, request.GET.get("speed_data").split(" "))
        full_results=""
        full_results = draw_all_by_speed(input_angles, speed_data, has_dispersia=False, full_results=full_results)
        fig = plt.gcf()
        fig.autofmt_xdate()
        plt.savefig(imgdata, format='png')
        imgdata.seek(0)
        print "Content-type: image/png\n"
        uri = 'data:image/png;base64,' + urllib.quote(base64.b64encode(imgdata.buf))
        imgdata.truncate()
        imgdata.seek(0)
        imgdata.close()
        # response = HttpResponse()
        # response.write(uri)
        response_data={
            "main_image": uri,
            "full_text": full_results
        }
        return HttpResponse(json.dumps(response_data), content_type="application/json")
