using System;

namespace ClassLibrary
{
    using System.Drawing;

    public class HatchingTransform : Interface.IPlugin
    {
        public string Name => "Штриховка";

        public string Version => "1.10";

        public string Author => "The best author";

        public void Transform(Interface.IMainApp app)
        {
            var bitmap = app.Image;
            for (var i = 0; i < bitmap.Width; i += 4)
                for (var j = 0; j < bitmap.Height; j++)
                    bitmap.SetPixel(i, j, Color.Black);

            app.Image = bitmap;
        }

    }
}
