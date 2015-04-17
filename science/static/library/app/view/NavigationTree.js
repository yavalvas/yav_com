Ext.define('Lib.view.NavigationTree', {
  extend: 'Ext.tree.Panel',
  alias: 'widget.navtree',
  initComponent: function() {
    var me = this,
      store = Ext.create('Ext.data.TreeStore', {
        with: 100,
        fields: [{
          name: 'text'
        }, {
          name: 'leaf'
        }, {
          name: 'tabId'
        }],
        root: {
          text: 'Линза Люнеберга ("Блитц")',
          expanded: true,
          children: [{
            text: "Модель линзы Люнеберга",
            leaf: true,
            tabId: 'lenseModel'
          }, {
            text: "График dl",
            leaf: true,
            tabId: 'dlGraph'
          },
          {
              text: "3D модель",
              leaf: true,
              tabId: "3dModel"
          },
          {
              text: "Описание модели",
              leaf: true,
              tabId: "descriptionModel"
          }
          ]
        }
      });

    Ext.apply(me, {
      store: store
    });

    me.callParent(arguments);
  }
});