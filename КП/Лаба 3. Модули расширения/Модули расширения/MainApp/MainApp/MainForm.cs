namespace MainApp
{
    using System;
    using System.Collections.Generic;
    using System.Drawing;
    using System.IO;
    using System.Reflection;
    using System.Windows.Forms;

    using Interface;

    public partial class MainForm : Form, IMainApp
    {
        private readonly Dictionary<string, IPlugin> plugins = new Dictionary<string, IPlugin>();

        public MainForm()
        {
            InitializeComponent();
            FindPlugins();
        }

        public Bitmap Image
        {
            get => (Bitmap)pcBox.Image;

            set => pcBox.Image = value;
        }

        private void FindPlugins()
        {
            var folder = AppDomain.CurrentDomain.BaseDirectory;
            var files = Directory.GetFiles(folder, "*.dll");
            foreach (var file in files)
            {
                try
                {
                    var assembly = Assembly.LoadFile(file);
                    foreach (var type in assembly.GetTypes())
                    {
                        var myInterface = type.GetInterface("Interface.IPlugin");
                        if (myInterface == null) continue;
                        var plugin = (IPlugin)Activator.CreateInstance(type);
                        plugins.Add(plugin.Name, plugin);
                    }
                }
                catch (Exception e)
                {
                    MessageBox.Show(@"Ошибка загрузки плагина" + e.Message);
                }
            }

            CreatePluginsMenu();
        }

        private void CreatePluginsMenu()
        {
            var item = new ToolStripMenuItem("Плагины");

            foreach (var name in plugins.Keys)
            {
                var pluginItem = new ToolStripMenuItem(name, null, OnPluginClick);
                item.DropDownItems.Add(pluginItem);
                filterMenu.Items.Add(item);
            }
        }

        private void OnPluginClick(object sender, EventArgs e)
        {
            var plugin = plugins[((ToolStripMenuItem)sender).Text];
            plugin.Transform(this);
        }
    }
}