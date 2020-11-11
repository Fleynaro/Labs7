using PluginInterface;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Imaging;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PluginContrast
{
    [Version(1, 0)]
    public class ContrastTransform : IPlugin
    {
        public string Name
        {
            get
            {
                return "Повысить контрастность";
            }
        }

        public string Author
        {
            get
            {
                return "Пухов Н. А.";
            }
        }

        public void Transform(Bitmap bitmap)
        {
            Bitmap originalImage = bitmap;
            Bitmap adjustedImage = (Bitmap)bitmap.Clone();
            float brightness = 1.0f; // no change in brightness
            float contrast = 1.5f; // twice the contrast
            float gamma = 1.0f; // no change in gamma

            float adjustedBrightness = brightness - 1.0f;
            // create matrix that will brighten and contrast the image
            float[][] ptsArray ={
                new float[] {contrast, 0, 0, 0, 0}, // scale red
                new float[] {0, contrast, 0, 0, 0}, // scale green
                new float[] {0, 0, contrast, 0, 0}, // scale blue
                new float[] {0, 0, 0, 1.0f, 0}, // don't scale alpha
                new float[] {adjustedBrightness, adjustedBrightness, adjustedBrightness, 0, 1}
            };

            ImageAttributes imageAttributes = new ImageAttributes();
            imageAttributes.ClearColorMatrix();
            imageAttributes.SetColorMatrix(new ColorMatrix(ptsArray), ColorMatrixFlag.Default, ColorAdjustType.Bitmap);
            imageAttributes.SetGamma(gamma, ColorAdjustType.Bitmap);
            Graphics g = Graphics.FromImage(adjustedImage); //холст
            g.DrawImage(originalImage, new Rectangle(0, 0, adjustedImage.Width, adjustedImage.Height)
                , 0, 0, originalImage.Width, originalImage.Height,
                GraphicsUnit.Pixel, imageAttributes); //рисуем на холсте, используя пиксели из originalImage с преобразованием imageAttributes

            for (int i = 0; i < bitmap.Width; ++i)
                for (int j = 0; j < bitmap.Height; ++j)
                {
                    Color color = adjustedImage.GetPixel(i, j);
                    bitmap.SetPixel(i, j, color);
                }
        }
    }

}
