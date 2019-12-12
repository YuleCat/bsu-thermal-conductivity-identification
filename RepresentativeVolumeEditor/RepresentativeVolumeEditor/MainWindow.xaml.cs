using RepresentativeVolumeEditor.Services;
using RepresentativeVolumeEditor.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace RepresentativeVolumeEditor
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        private List<Material> materials;
        private Label activeColorLabel;

        public MainWindow()
        {
            InitializeComponent();

            this.materials = RepresentativeVolumeMaterialsService.GetMaterials();
            InitializeMaterialsGrid();

            RepaintRepresentativeVolume();
        }

        private void InitializeMaterialsGrid()
        {
            var materialsCount = this.materials.Count;
            for (int i = 0; i < materialsCount; ++i)
            {
                this.materialsGrid.RowDefinitions.Add(new RowDefinition());

                var colorLabel = CreateColorLabel(this.materials[i].Color);
                this.materialsGrid.Children.Add(colorLabel);
                Grid.SetRow(colorLabel, i);
                Grid.SetColumn(colorLabel, 0);

                var nameLabel = CreateNameLabel(this.materials[i].Name);
                this.materialsGrid.Children.Add(nameLabel);
                Grid.SetRow(nameLabel, i);
                Grid.SetColumn(nameLabel, 1);
            }
        }

        private void PeriodTextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            var rawPeriod = this.periodTextBox.Text;
            this.periodTextBox.BorderThickness = new Thickness(1);
            if (!int.TryParse(rawPeriod, out int period) || period <= 0 || period >= 500)
            {
                this.periodTextBox.BorderBrush = new SolidColorBrush(Colors.Red);
                return;
            }

            this.periodTextBox.BorderBrush = new SolidColorBrush(Colors.Black);

            RepaintRepresentativeVolume();
        }

        private void CellSizeTextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            var rawCellSize = this.cellSizeTextBox.Text;
            this.cellSizeTextBox.BorderThickness = new Thickness(1);
            if (!int.TryParse(rawCellSize, out int cellSize) || cellSize < 5 || cellSize > 250)
            {
                this.cellSizeTextBox.BorderBrush = new SolidColorBrush(Colors.Red);
                return;
            }

            this.cellSizeTextBox.BorderBrush = new SolidColorBrush(Colors.Black);

            RepaintRepresentativeVolume();
        }

        private void ColorLabel_MouseLeftButtonUp(object sender, MouseButtonEventArgs e)
        {
            var colorLabel = (Label)sender;
            if (colorLabel == activeColorLabel)
            {
                return;
            }
            colorLabel.BorderBrush = new SolidColorBrush(Colors.Red);
            if (activeColorLabel != null)
            {
                activeColorLabel.BorderBrush = new SolidColorBrush(Colors.Black);
            }
            activeColorLabel = colorLabel;
        }

        private void ColorLabel_MouseLeave(object sender, MouseEventArgs e)
        {
            if (this.Cursor != Cursors.Wait)
            {
                Mouse.OverrideCursor = Cursors.Arrow;
            }
        }

        private void ColorLabel_MouseEnter(object sender, MouseEventArgs e)
        {
            if (this.Cursor != Cursors.Wait)
            {
                Mouse.OverrideCursor = Cursors.Hand;
            }
        }

        private void Cell_MouseLeftButtonUp(object sender, MouseButtonEventArgs e)
        {
            if (activeColorLabel == null)
            {
                return;
            }

            var cell = (Label)sender;
            if ((cell.Background as SolidColorBrush).Color == (activeColorLabel.Background as SolidColorBrush).Color)
            {
                cell.Background = new SolidColorBrush(Colors.Gray);
                return;
            }
            cell.Background = activeColorLabel.Background;
        }

        private void SaveButton_MouseLeftButtonUp(object sender, EventArgs e)
        {
            var grid = new List<List<Label>>();
            for (int i = 0; i < volumeGrid.Children.Count; ++i)
            {
                grid.Add(new List<Label>());
                var row = volumeGrid.Children[i] as Grid;
                for (int j = 0; j < row.Children.Count; ++j)
                {
                    var label = row.Children[j] as Label;
                    if (label != null)
                    {
                        grid[i].Add(label);
                    }
                }
            }

            var period = int.Parse(this.periodTextBox.Text);
            var cellSize = int.Parse(this.cellSizeTextBox.Text);

            RepresentativeVolumeSaverService.Save(materials, grid, period, cellSize, (bool)this.chessOrderCheckBox.IsChecked);
        }

        private void ChessOrderCheckBox_Click(object sender, RoutedEventArgs e)
        {
            RepaintRepresentativeVolume();
        }

        private void RepaintRepresentativeVolume()
        {
            if (this.periodTextBox == null || this.cellSizeTextBox == null || this.chessOrderCheckBox == null)
            {
                return;
            }

            var period = int.Parse(this.periodTextBox.Text);
            var cellSize = int.Parse(this.cellSizeTextBox.Text);

            if (cellSize > period)
            {
                this.periodTextBox.BorderBrush = new SolidColorBrush(Colors.Red);
                this.cellSizeTextBox.BorderBrush = new SolidColorBrush(Colors.Red);
                return;
            }

            this.periodTextBox.BorderBrush = new SolidColorBrush(Colors.Black);
            this.cellSizeTextBox.BorderBrush = new SolidColorBrush(Colors.Black);

            this.volumeGrid.Children.Clear();
            this.volumeGrid.RowDefinitions.Clear();

            var isChessOrder = (bool)this.chessOrderCheckBox.IsChecked;
            var oddCellCount = (int)Math.Ceiling(this.volumeGrid.Width / period);
            var evenCellCount = (int)Math.Ceiling((this.volumeGrid.Width - period / 2) / period);

            for (int i = 0; i < oddCellCount; ++i)
            {
                var row = new RowDefinition();
                row.Height = new GridLength(period);
                this.volumeGrid.RowDefinitions.Add(row);

                var rowGrid = new Grid();
                var innerRow = new RowDefinition();
                rowGrid.RowDefinitions.Add(innerRow);
                if (!isChessOrder || isChessOrder && (i + 1) % 2 != 0)
                {
                    for (int j = 0; j < oddCellCount; ++j)
                    {
                        var column = new ColumnDefinition();
                        column.Width = new GridLength(period);
                        rowGrid.ColumnDefinitions.Add(column);

                        var cell = CreateCell(cellSize);
                        rowGrid.Children.Add(cell);
                        Grid.SetRow(cell, 0);
                        Grid.SetColumn(cell, j);
                    }
                }
                else
                {
                    var extraColumn = new ColumnDefinition();
                    extraColumn.Width = new GridLength(period / 2);
                    rowGrid.ColumnDefinitions.Add(extraColumn);

                    for (int j = 1; j <= evenCellCount; ++j)
                    {
                        var column = new ColumnDefinition();
                        column.Width = new GridLength(period);
                        rowGrid.ColumnDefinitions.Add(column);

                        var cell = CreateCell(cellSize);
                        rowGrid.Children.Add(cell);
                        Grid.SetRow(cell, 0);
                        Grid.SetColumn(cell, j);
                    }
                }
                volumeGrid.Children.Add(rowGrid);
                Grid.SetRow(rowGrid, i);
            }
        }

        private Label CreateCell(int cellSize)
        {
            var cell = new Label();
            cell.Background = new SolidColorBrush(Colors.Gray);
            cell.BorderBrush = new SolidColorBrush(Colors.Black);
            cell.BorderThickness = new Thickness(1);
            cell.Width = cellSize;
            cell.Height = cellSize;
            cell.HorizontalAlignment = HorizontalAlignment.Left;
            cell.VerticalAlignment = VerticalAlignment.Top;
            cell.MouseEnter += ColorLabel_MouseEnter;
            cell.MouseLeave += ColorLabel_MouseLeave;
            cell.MouseLeftButtonUp += Cell_MouseLeftButtonUp;

            return cell;
        }

        private Label CreateColorLabel(Color color)
        {
            var colorLabel = new Label();
            colorLabel.Background = new SolidColorBrush(color);
            colorLabel.BorderBrush = new SolidColorBrush(Colors.Black);
            colorLabel.BorderThickness = new Thickness(2);
            colorLabel.Width = 20;
            colorLabel.Height = 20;
            colorLabel.Margin = new Thickness(5, 0, 0, 2);
            colorLabel.MouseEnter += ColorLabel_MouseEnter;
            colorLabel.MouseLeave += ColorLabel_MouseLeave;
            colorLabel.MouseLeftButtonUp += ColorLabel_MouseLeftButtonUp;

            return colorLabel;
        }

        private Label CreateNameLabel(string materialName)
        {
            var nameLabel = new Label();
            nameLabel.Content = materialName;

            return nameLabel;
        }
    }
}
