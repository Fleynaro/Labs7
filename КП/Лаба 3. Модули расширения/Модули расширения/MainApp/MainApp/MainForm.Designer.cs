namespace MainApp
{
    partial class MainForm
    {
        /// <summary>
        /// Обязательная переменная конструктора.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Освободить все используемые ресурсы.
        /// </summary>
        /// <param name="disposing">истинно, если управляемый ресурс должен быть удален; иначе ложно.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Код, автоматически созданный конструктором форм Windows

        /// <summary>
        /// Требуемый метод для поддержки конструктора — не изменяйте 
        /// содержимое этого метода с помощью редактора кода.
        /// </summary>
        private void InitializeComponent()
        {
            this.pcBox = new System.Windows.Forms.PictureBox();
            this.filterMenu = new System.Windows.Forms.MenuStrip();
            this.openImageDialog = new System.Windows.Forms.OpenFileDialog();
            this.saveImageDialog = new System.Windows.Forms.SaveFileDialog();
            ((System.ComponentModel.ISupportInitialize)(this.pcBox)).BeginInit();
            this.SuspendLayout();
            // 
            // pcBox
            // 
            this.pcBox.BackgroundImage = global::MainApp.Properties.Resources.MainImage;
            this.pcBox.BackgroundImageLayout = System.Windows.Forms.ImageLayout.Center;
            this.pcBox.Dock = System.Windows.Forms.DockStyle.Fill;
            this.pcBox.Image = global::MainApp.Properties.Resources.MainImage;
            this.pcBox.InitialImage = null;
            this.pcBox.Location = new System.Drawing.Point(0, 24);
            this.pcBox.Name = "pcBox";
            this.pcBox.Size = new System.Drawing.Size(786, 490);
            this.pcBox.SizeMode = System.Windows.Forms.PictureBoxSizeMode.StretchImage;
            this.pcBox.TabIndex = 0;
            this.pcBox.TabStop = false;
            // 
            // filterMenu
            // 
            this.filterMenu.Location = new System.Drawing.Point(0, 0);
            this.filterMenu.Name = "filterMenu";
            this.filterMenu.Size = new System.Drawing.Size(786, 24);
            this.filterMenu.TabIndex = 1;
            this.filterMenu.Text = "menuStrip1";
            // 
            // openImageDialog
            // 
            this.openImageDialog.FileName = "openFileDialog1";
            // 
            // MainForm
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(786, 514);
            this.Controls.Add(this.pcBox);
            this.Controls.Add(this.filterMenu);
            this.FormBorderStyle = System.Windows.Forms.FormBorderStyle.FixedSingle;
            this.MainMenuStrip = this.filterMenu;
            this.MaximizeBox = false;
            this.Name = "MainForm";
            this.Text = "Фильтры";
            ((System.ComponentModel.ISupportInitialize)(this.pcBox)).EndInit();
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.PictureBox pcBox;
        private System.Windows.Forms.MenuStrip filterMenu;
        private System.Windows.Forms.OpenFileDialog openImageDialog;
        private System.Windows.Forms.SaveFileDialog saveImageDialog;
    }
}

