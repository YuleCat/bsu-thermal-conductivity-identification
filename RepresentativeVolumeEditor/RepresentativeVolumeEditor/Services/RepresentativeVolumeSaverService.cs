using RepresentativeVolumeEditor.Utils;
using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using System.Windows.Controls;
using System.Windows.Media;

namespace RepresentativeVolumeEditor.Services
{
    static class RepresentativeVolumeSaverService
    {
        public static void Save(List<Material> materials, List<List<Label>> grid, int period, int cellSize, bool isChessOrder)
        {
            var materialsMappingBuilder = new StringBuilder();
            var colorsMapping = new Dictionary<Color, int>();
            for (int i = 0; i < materials.Count; ++i)
            {
                materialsMappingBuilder.Append($"{materials[i].Name}={i};");
                colorsMapping.Add(materials[i].Color, i);
            }
            colorsMapping.Add(Colors.Gray, colorsMapping.Count);

            var savePath = "./data/representative_volume.txt";
            Directory.CreateDirectory(Path.GetDirectoryName(savePath));
            using (var file = new StreamWriter(savePath))
            {
                file.WriteLine($"period={period};");
                file.WriteLine($"cell_size={cellSize};");
                file.WriteLine($"chess_order={isChessOrder.ToString()};");
                file.WriteLine($"materials:{materialsMappingBuilder.ToString()}");
                var rowInfoBuilder = new StringBuilder();
                foreach (var row in grid)
                {
                    foreach (var cell in row)
                    {
                        var cellColor = (cell.Background as SolidColorBrush).Color;
                        rowInfoBuilder.Append($"{colorsMapping[cellColor]} ");
                    }
                    rowInfoBuilder.Append("\n");
                }
                file.WriteLine($"grid:");
                file.WriteLine(rowInfoBuilder.ToString());
            }
        }
    }
}
