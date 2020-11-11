namespace ClassLibrary
{
    using System.Drawing;

    using Interface;

    public class HalfCoupTransform : IPlugin
    {
        public string Name => "Замена красного и синего цветов";

        public string Version => "1.10";

        public string Author => "The best author";

        public void Transform(IMainApp app)
        {
            var bitmap = app.Image;
            for (var i = 0; i < bitmap.Width; i++)
            {
                for (var j = 0; j < bitmap.Height; j++)
                {
                    var currentColor = bitmap.GetPixel(i, j);
                    bitmap.SetPixel(i, j, Color.FromArgb(currentColor.B, currentColor.G, currentColor.R));
                }
            }

            app.Image = bitmap;
        }
    }
}