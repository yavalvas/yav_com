function detectmobBySize() {
    if (window.innerWidth <= 800 && window.innerHeight <= 600) {
        return true;
    } else {
        return false;
    }
}
function detectmob() {
    if (navigator.userAgent.match(/Android/i)
        || navigator.userAgent.match(/webOS/i)
        || navigator.userAgent.match(/iPhone/i)
        || navigator.userAgent.match(/iPad/i)
        || navigator.userAgent.match(/iPod/i)
        || navigator.userAgent.match(/BlackBerry/i)
        || navigator.userAgent.match(/Windows Phone/i)
    ) {
        return true;
    }
    else {
        return false;
    }
}
Ext.tip.QuickTipManager.init();
Ext.application({
    name: 'Lib',
    appFolder: '/static/library/app',
    controllers: ['Main'],
    views: ['Main'],
    //views: ['Main', 'Ext.ux.panel.PDF'],
//  models: ['Genre'],

    launch: function () {
        var mainView = Ext.create('Ext.container.Viewport', {
            layout: 'fit',
            items: {
                xtype: 'mainview'
            }

        });
        var isMobileBySize = detectmobBySize();
        console.log(isMobileBySize);
        if (isMobileBySize) {
            Ext.getCmp("mainTree").collapse();
        }
        var isMobile = detectmob();
        console.log(isMobile);
        if (isMobile) {
            Ext.getCmp("mainTree").collapse();
        }

        Ext.EventManager.onWindowResize(function (w, h) {
            console.log(1);
            isMobileBySize = detectmobBySize();
            if (detectmobBySize()) {
                Ext.getCmp("mainTree").collapse();
            }
            else {
                Ext.getCmp("mainTree").expand();
            }
            //mainView.doComponentLayout();

        });


    }
});