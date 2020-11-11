using PluginInterface;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.IO;
using System.Reflection;
using System.Windows.Forms;
using System.Xml;

namespace ExtModule
{
    public partial class Form1 : Form
    {
        List<string> foundPlugins = new List<string>();
       
        struct PluginInfo
        {
            public IPlugin plugin;
            public VersionAttribute version;

            public PluginInfo(IPlugin plugin, VersionAttribute version)
            {
                this.plugin = plugin;
                this.version = version;
            }
        }

        Dictionary<string, PluginInfo> plugins = new Dictionary<string, PluginInfo>();
       
        public Form1()
        {
            InitializeComponent();
            FindPlugins();
            LoadFoundPlugins();
            CreatePluginsMenu();
        }

        private void LoadFoundPlugins()
        {
            foreach (string file in foundPlugins)
            {
                Assembly assembly = Assembly.LoadFrom(file);

                foreach (Type type in assembly.GetTypes())
                {
                    Type iface = type.GetInterface("PluginInterface.IPlugin");

                    if (iface != null)
                    {
                        var obj = Activator.CreateInstance(type);
                        var plugin = (IPlugin)obj;
                        var attr = (VersionAttribute)type.GetCustomAttribute(typeof(VersionAttribute));
                        plugins.Add(plugin.Name, new PluginInfo(plugin, attr));                   
                    }
                }
            }
        }

        private void FindPlugins()
        {
            // папка с плагинами
            string folder = AppDomain.CurrentDomain.BaseDirectory;
            var domain = AppDomain.CreateDomain("temp");
            try
            {
                XmlDocument xDoc = new XmlDocument();
                xDoc.Load(folder + "plugins.xml");

                XmlElement tagPlugins = xDoc.DocumentElement;

                //load from cfg file
                XmlNode attr = tagPlugins.Attributes.GetNamedItem("loading");
                if (attr.Value == "manual")
                {
                    foreach (XmlNode tagPlugin in tagPlugins)
                    {
                        FindPlugin(folder + tagPlugin.Attributes.GetNamedItem("path").Value, domain);
                    }
                    AppDomain.Unload(domain);
                    return;
                }
            } catch(System.IO.FileNotFoundException) { }

            //load all
            foreach (string file in Directory.GetFiles(folder, "*.dll"))
            {
                FindPlugin(file, domain);
            }
            AppDomain.Unload(domain);
        }

        private void FindPlugin(string file, AppDomain domain)
        {
            try
            {
                Assembly assembly = domain.Load(File.ReadAllBytes(file));

                foreach (Type type in assembly.GetTypes())
                {
                    Type iface = type.GetInterface("PluginInterface.IPlugin");

                    if (iface != null)
                    {
                        foundPlugins.Add(file);
                    }
                }
            }
            catch (Exception ex)
            {
                MessageBox.Show("Ошибка загрузки плагина\n" + ex.Message);
            }
        }

        private void CreatePluginsMenu()
        {
            var item = new ToolStripMenuItem("Filters");
            menuStrip1.Items.Add(item);

            foreach (var name in plugins.Keys)
            {
                var pluginItem = new ToolStripMenuItem(name, null, OnPluginClick);
                item.DropDownItems.Add(pluginItem);
            }

            item.DropDownItems.Add(new ToolStripMenuItem("Список", null, OnPluginListClick));
        }

        private void OnPluginClick(object sender, EventArgs e)
        {
            IPlugin plugin = plugins[((ToolStripMenuItem)sender).Text].plugin;
            plugin.Transform((Bitmap)pictureBox1.Image);
            pictureBox1.Invalidate();
        }

        private void OnPluginListClick(object sender, EventArgs e)
        {
            string info = "";
            foreach (var plugin in plugins)
            {
                info += plugin.Value.plugin.Name + " (Author: "+ plugin.Value.plugin.Author + " | Version: "+ plugin.Value.version.Major + "." + plugin.Value.version.Minor + ")\n";
            }
            MessageBox.Show(info);
        }
    }
}
