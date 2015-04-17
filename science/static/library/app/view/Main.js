function base64toBlob(base64Data, contentType) {
    contentType = contentType || '';
    var sliceSize = 1024;
    var byteCharacters = atob(base64Data);
    var bytesLength = byteCharacters.length;
    var slicesCount = Math.ceil(bytesLength / sliceSize);
    var byteArrays = new Array(slicesCount);

    for (var sliceIndex = 0; sliceIndex < slicesCount; ++sliceIndex) {
        var begin = sliceIndex * sliceSize;
        var end = Math.min(begin + sliceSize, bytesLength);

        var bytes = new Array(end - begin);
        for (var offset = begin, i = 0; offset < end; ++i, ++offset) {
            bytes[i] = byteCharacters[offset].charCodeAt(0);
        }
        byteArrays[sliceIndex] = new Uint8Array(bytes);
    }
    return new Blob(byteArrays, {type: contentType});
}
//???
//var pathToImg = {% MEDIA animated.gif %}

var pdfPanel = Ext.create('Ext.ux.panel.PDF', {
        title    : 'PDF Panel',
        width    : 489,
        height   : 633,
        pageScale: 0.75,                                           // Initial scaling of the PDF. 1 = 100%
        src      : 'media/tracemonkey.pdf' // URL to the PDF - Same Domain or Server with CORS Support
        //renderTo : Ext.get('pdfPanel')
    });
Ext.define('Lib.view.Main', {
    extend: 'Ext.container.Container',
    frame: true,
    border: true,
    closable: true,
    flex: 1,
    region: 'center',
    width: "100%",
    height: "100%",
    requires: [
        'Ext.grid.Panel',
        'Ext.tab.Panel',
        'Ext.button.Button',
        'Lib.view.NavigationTree',
        "Ext.ux.panel.PDF"
    ],
    alias: 'widget.mainview',
    layout: {
        type: 'border',
        align: 'stretch'
    },
    items: [
        {
            xtype: 'navtree',
            region: 'west',
            width: 190,
            id: "mainTree",
            collapsible: true,
            collapseMode: 'mini'
        },
        {
            xtype: 'tabpanel',
            region: 'center',
            items: [
                {
                    title: 'Модель линзы',
                    xtype: 'panel',
                    id: "lenseModel",
                    layout: {
                        type: 'hbox',
                        align: 'stretch'
                    },
                    items: [
                        {
                            xtype: 'form',
                            collapsible: true,
                            id: 'idForm',
                            frame: true,
                            width: 200,
                            height: 300,
                            //layout: {
                            //    type: 'absolute'
                            //},
                            bodyPadding: 10,
                            title: 'Расчеты хода луча',
                            items: [
                                {
                                    xtype: 'numberfield',
                                    id: 'idVal1',
                                    fieldLabel: 'Vo [м/с]',
                                    value: 7500,
                                    width: 180,
                                    x: 10,
                                    y: 10
                                },
                                {
                                    xtype: 'numberfield',
                                    id: 'idVal2',
                                    fieldLabel: 'R1 [м]',
                                    width: 180,
                                    value: 0.085,
                                    x: 10,
                                    y: 35,
                                    decimalPrecision: 5
                                },
                                {
                                    xtype: 'numberfield',
                                    id: 'idVal3',
                                    fieldLabel: 'R2 [м]',
                                    width: 180,
                                    value: 0.0535,
                                    x: 10,
                                    y: 60,
                                    decimalPrecision: 5
                                },
                                {
                                    xtype: 'numberfield',
                                    id: 'idVal4',
                                    fieldLabel: 'n1 [ед]',
                                    width: 180,
                                    value: 1.47290,
                                    x: 10,
                                    y: 85,
                                    decimalPrecision: 5
                                },
                                {
                                    xtype: 'numberfield',
                                    id: 'idVal5',
                                    fieldLabel: 'n2 [ед]',
                                    width: 180,
                                    value: 1.76470,
                                    x: 10,
                                    y: 110,
                                    decimalPrecision: 5
                                },
                                {
                                    xtype: 'numberfield',
                                    id: 'idVal6',
                                    width: 180,
                                    fieldLabel: 'R_out[м]',
                                    value: 0.12,
                                    x: 10,
                                    y: 135,
                                    decimalPrecision: 5
                                },
                                {
                                    xtype: 'textfield',
                                    id: 'idAngles',
                                    fieldLabel: 'Углы [град]',
                                    defaultValue: "215 220 225 230",
                                    width: 180,
                                    value: "215 220 225 230",
                                    x: 10,
                                    y: 160
                                },
                                {
                                    xtype: "checkbox",
                                    x: 10,
                                    y: 185,
                                    width: 180,
                                    id: "withSharp",
                                    fieldLabel: "Показать сетку"
                                },
                                {
                                    xtype: "checkbox",
                                    x: 10,
                                    y: 210,
                                    width: 180,
                                    id: "withCoords",
                                    fieldLabel: "Показать координаты точек"
                                },
                                {
                                    xtype: "checkbox",
                                    x: 10,
                                    y: 245,
                                    width: 180,
                                    id: "withRaws",
                                    fieldLabel: "Показать со стрелками"
                                },
                                {
                                    xtype: "checkbox",
                                    x: 10,
                                    y: 270,
                                    id: "hasDispersia",
                                    width: 180,
                                    fieldLabel: "Расчет с дисперсией"
                                },
                                {
                                    xtype: 'button',
                                    text: 'Рассчитать',
                                    region: "north",
                                    handler: function () {
//                              var lVal1 = Ext.getCmp('idVal1').value,
//                                  lVal2 = Ext.getCmp('idVal2').value,
//                                  lOp = Ext.getCmp('idOp'),
//                                  lResult = Ext.getCmp('idResult'),
//                                  lRecord = lOp.findRecordByDisplay(lOp.rawValue);

//                              DemoEasyExtJS4.Compute.Execute(lVal1,lRecord.raw.exec,lVal2, function(pResult){
//                                lResult.setValue(pResult);
//                              });
                                        var Mask = new Ext.LoadMask(Ext.getCmp("displayImg"), {msg: "Обновление построения"});
                                        Mask.show();
                                        Ext.Ajax.request({
                                            url: '/graphics/graphic',
                                            method: "GET",
                                            params: {
                                                Vo: Ext.getCmp("idVal1").getValue().toString(),
                                                R1: Ext.getCmp("idVal2").getValue().toString(),
                                                R2: Ext.getCmp("idVal3").getValue().toString(),
                                                N1: Ext.getCmp("idVal4").getValue().toString(),
                                                N2: Ext.getCmp("idVal5").getValue().toString(),
                                                R_out: Ext.getCmp("idVal6").getValue().toString(),
                                                angles: Ext.getCmp("idAngles").getValue().toString(),
                                                withSharp: Ext.getCmp("withSharp"),
                                                withCoords: Ext.getCmp("withCoords"),
                                                withArrows: Ext.getCmp("withRaws"),
                                                hasDispersia: Ext.getCmp("hasDispersia").getValue().toString()
                                            },
                                            success: function (response, opts) {
                                                var full_info = JSON.parse(response.responseText);
                                                if(full_info.hasOwnProperty("response_error")){
                                                    switch(full_info.response_error){
                                                        case "wrong_angle_type":
                                                            Ext.Msg.alert("Сообщение", "Неверно задан угол");
                                                            break;
                                                        case "wrong_angle_value":
                                                            Ext.Msg.alert("Сообщение", "Значения углов должны находиться в пределах от 215 до 325");
                                                            break;
                                                        case "r1_lower_r2":
                                                            Ext.Msg.alert("Сообщение", "Радиус R1 меньше, чем R2");
                                                            break;
                                                        case "wrong_vo":
                                                            Ext.Msg.alert("Сообщение", "Неверное значение скорости");
                                                            break;
                                                        case "wrong_value":
                                                            Ext.Msg.alert("Сообщение", "Неверные значения среди R1, R2, N1, N2");
                                                            break;
                                                    }
                                                }
                                                else {
                                                    Ext.getCmp("displayImg").setSrc(full_info.main_image);
                                                    Ext.getCmp("displayOutput").setValue(full_info.full_text);
                                                }
                                                Mask.hide();
                                            },
                                            failure: function (form, action) {
                                                console.log("Произошла ошибка");
                                                Ext.Msg.alert("Сообщение", "Произошла ошибка");
                                                Mask.hide();
                                            }
                                        })
                                    }
                                    //x: 180,
                                    //y: 250
                                }
                            ]
                        },
                        {
                            xtype: 'panel',
                            title: 'Построение',
                            width: "auto",
                            height: "auto",
                            id: 'graphPnl',
                            region: "north",
                            collapsible: true,
                            autoScroll: true,
                            items: [
                                {
                                    xtype: 'image',
                                    id: 'displayImg',
                                    //width: 600,
                                    width: 500,
                                    //height: 600,
                                    height: 500,
                                    style: {
                                        'display': 'block',
                                        'margin': 'auto'
                                    },
//                            renderTo: Ext.getBody(),
//              src: "/graphics/graphic.png",
                                    region: "north"
                                },
                                {
                                    xtype: "displayfield",
                                    id: "displayOutput",
                                    width: 600,
                                    height: 300,
                                    region: "south",
                                    renderer: function (value, x) {
                                        //hack for displayfield issue in ExtJS 5
                                        var rtn = value.replace(/(\r\n|\n|\r)/gm, "</br>");
                                        return rtn
                                    }
                                }
                            ]
                        }
                    ]
                },
                {
                    title: 'График dl',
                    itemId: 'dlGraph',
                    xtype: "panel",
                    layout: {
                        type: 'hbox',
                        align: 'stretch'
                    },
                    items: [
                        {
                            xtype: 'form',
                            collapsible: true,
                            id: 'idFormDl',
                            frame: true,
                            width: 270,
                            height: 300,
                            bodyPadding: 10,
                            title: 'Расчеты отклонения dl',
                            items: [
                                {
                                    xtype: "checkbox",
                                    x: 10,
                                    y: 10,
                                    id: "withSharpDl",
                                    fieldLabel: "Показать сетку"
                                },
                                {
                                    xtype: 'textfield',
                                    id: 'idSpeeds',
                                    fieldLabel: 'Скорости Vo [м/с]',
                                    defaultValue: "0 7500 10000 20000",
                                    value: "0 7500 10000 20000",
                                    x: 10,
                                    y: 35
                                },
                                {
                                    xtype: 'textfield',
                                    id: 'idAnglesDl',
                                    fieldLabel: 'Углы [град]',
                                    defaultValue: "225 230 235 240 245 250 255 260 265 267 268 269 271 272 273 274 275 280 285 290 295 300 310 315",
                                    value: "225 230 235 240 245 250 255 260 265 267 268 269 271 272 273 274 275 280 285 290 295 300 310 315",
                                    x: 10,
                                    y: 60
                                },
                                {
                                    xtype: 'button',
                                    id: 'idPlotDiff',
                                    text: 'Построить',
                                    x: 10,
                                    y: 85,
                                    handler: function () {
                                        var Mask = new Ext.LoadMask(Ext.getCmp("displayImgDl"), {msg: "Обновление построения"});
                                        Mask.show();
                                        Ext.Ajax.request({
                                            url: '/graphics/plot_diff',
                                            method: "GET",
                                            params: {
                                                withSharp: Ext.getCmp("withSharpDl")
                                            },
                                            success: function (response, opts) {
                                                var full_info = JSON.parse(response.responseText);
                                                console.log(full_info);
                                                Ext.getCmp("displayImgDl").setSrc(full_info.main_image);
                                                Ext.getCmp("displayOutputDl").setValue(full_info.full_text);
                                                Mask.hide();
                                            },
                                            failure: function (form, action) {
                                                console.log("Произошла ошибка");
                                                Ext.Msg.alert("Сообщение", "Произошла ошибка");
                                                Mask.hide();
                                            }
                                        })
                                    }
                                },
                                {
                                    xtype: 'button',
                                    id: 'idOutputDiff',
                                    text: 'Построить разность',
                                    x: 10,
                                    y: 120,
                                    handler: function () {
                                        var Mask = new Ext.LoadMask(Ext.getCmp("displayImgDl"), {msg: "Обновление построения"});
                                        Mask.show();
                                        Ext.Ajax.request({
                                            url: '/graphics/all_by_speed',
                                            method: "GET",
                                            params: {
                                                withSharp: Ext.getCmp("withSharpDl"),
                                                angles: Ext.getCmp("idAnglesDl").getValue().toString(),
                                                speed_data: Ext.getCmp("idSpeeds").getValue().toString()
                                            },
                                            success: function (response, opts) {
                                                var full_info = JSON.parse(response.responseText);
                                                console.log(full_info);
                                                Ext.getCmp("displayImgDl").setSrc(full_info.main_image);
                                                Ext.getCmp("displayOutputDl").setValue(full_info.full_text);
                                                Mask.hide();
                                            },
                                            failure: function (form, action) {
                                                console.log("Произошла ошибка");
                                                Ext.Msg.alert("Сообщение", "Произошла ошибка");
                                                Mask.hide();
                                            }
                                        })
                                    }
                                }
                            ]
                        },
                        {
                            xtype: 'panel',
                            title: 'Построение',
                            width: "auto",
                            height: "auto",
                            id: 'graphDlPnl',
                            region: "north",
                            collapsible: true,
                            autoScroll: true,
                            items: [
                                {
                                    xtype: 'image',
                                    id: 'displayImgDl',
                                    //width: 600,
                                    width: 500,
                                    //height: 600,
                                    height: 500,
                                    style: {
                                        'display': 'block',
                                        'margin': 'auto'
                                    },
//                            renderTo: Ext.getBody(),
//              src: "/graphics/graphic.png",
                                    region: "north"
                                },
                                {
                                    xtype: "displayfield",
                                    id: "displayOutputDl",
                                    width: 600,
                                    height: 300,
                                    region: "south",
                                    renderer: function (value, x) {
                                        //hack for displayfield issue in ExtJS 5
                                        var rtn = value.replace(/(\r\n|\n|\r)/gm, "</br>");
                                        return rtn
                                    }
                                }
                            ]
                        }
                    ]
                },
                {
                    title: 'Трехмерная модель',
                    itemId: '3dModel',
                    xtype: 'panel',
                    items: [
                        {
                            xtype: "numberfield",
                            minValue: 0,
                            maxValue: 299792458,
                            width: 45,
                            value: 7500
                        },
                        {
                            xtype: "slider",
                            fieldLabel: "Скорость",
                            minValue: 0,
                            maxValue: 299792458,
                            width: 450,
                            value: 7500
                        },
                        {
                            xtype: "panel",
                            id: '3dviewer',
                            html: '<canvas id="canvas" style="border:5px solid #000 width="700" height="700"></canvas>',
                            //autoEl: {
                            //    tag: 'canvas',
                            //    height: 600,
                            //    width: 600
                            //},
                            width: 700,
                            height: 700,
                            //width: "100%",
                            //height: "100%",
                            listeners: {
                                'afterlayout': function (form) {
                                    var camera, scene, material, mesh, geometry, renderer;

                                    function drawSphere() {
                                        init();
                                        animate();

                                    }

                                    function init() {
                                        // camera
                                        //container = document.getElementById("containerSphere");
                                        scene = new THREE.Scene();
                                        canvas = document.getElementById("canvas");
                                        camera = new THREE.PerspectiveCamera(50, window.innerWidth / innerHeight, 1, 1000);
                                        camera.position.z = 300;
                                        scene.add(camera);

                                        // sphere object
                                        var radius = 50,
                                            segments = 10,
                                            rings = 10;
                                        geometry = new THREE.SphereGeometry(radius, segments, rings);
                                        material = new THREE.MeshNormalMaterial({color: 0x002288});
                                        mesh = new THREE.Mesh(geometry, material);
                                        //scene
                                        scene.add(mesh);


                                        //renderer.setSize(300, 300);
                                        //Ext.getCmp('3dviewer').body.appendChild(renderer.domElement);
                                        //Ext.getCmp('3dviewer').add(Ext.create('Ext.Panel', {
                                        //    html: '<div id="canvasContainer"></div>'
                                        //}));
                                        //Ext.getCmp('3dviewer').add(Ext.create('Ext.Panel', {
                                        //    id: "canvasElement",
                                        //    html: '<canvas id="myCanvas"></canvas>'
                                        //}));

                                        // renderer
                                        //var canvas = Ext.getCmp("3dviewer").html;
                                        //renderer = new THREE.WebGLRenderer({
                                        //antialias: true
                                        //canvas: canvas
                                        //});
                                        co = canvas.getContext("2d");
                                        renderer = new THREE.CanvasRenderer({canvas: canvas});

                                        renderer.setSize(window.innerWidth, window.innerHeight);        //Ext.getCmp('3dviewer').html="'"+renderer.domElement+"'";
                                        //Ext.getCmp('3dviewer').add(Ext.create('Ext.Panel', {
                                        //    html: renderer.domElement
                                        //}));
                                        //document.body.appendChild(renderer.domElement);
                                        //var c = document.getElementById("myCanvas");
                                        //renderer = new THREE.CanvasRenderer({canvas: c});
                                        //renderer.setSize(700, 700);

                                    }


                                    function animate() {
                                        requestAnimationFrame(animate);
                                        co.clearRect(0, 0, canvas.width, canvas.height);
                                        render();

                                    }

                                    function render() {

                                        mesh.rotation.x += .01;
                                        mesh.rotation.y += .02;
                                        renderer.render(scene, camera);


                                    }

                                    // fn callin
                                    drawSphere();
                                }
                            }
                            //region: 'center',
                            //layout: 'fit'

                        }
                        //{
                        //    xtype: 'box',
                        //    autoEl:{
                        //        tag: 'canvas',
                        //        height: 600,
                        //        width: 600,
                        //        id: "canvasElement"
                        //    }
                        //}
                        //{xtype: "component",
                        // width: 576,
                        // height: 420,
                        // autoEl: {
                        //    tag : "iframe"
                        //    ,src: "http://slides.com/yavalvas/china_economy/embed"
                        //    ,frameborder: 0,
                        //     width: 576,
                        //     height: 420,
                        //     scrolling: "no",
                        //     frameborder: "0",
                        //     webkitallowfullscreen: true,
                        //     mozallowfullscreen: true,
                        //     allowfullscreen: true
                        //}}
                    ]
                },
                {
                    title: "Описание расчетов",
                    itemId: 'descriptionModel',
                    xtype: 'panel',
                    items: [
                        pdfPanel
                    ]
                }
            ]
        }
    ]
});