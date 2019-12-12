using MongoDB.Bson;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using System.Windows.Media;

namespace RepresentativeVolumeEditor.Utils
{
    class Material
    {
        public Material(string name, Color color, double thermalConductivity, double density, double heatCapacity)
        {
            if (string.IsNullOrEmpty(name))
            {
                throw new ArgumentNullException(nameof(name));
            }

            this.Name = name;
            this.Color = color;
            this.ThermalConductivity = thermalConductivity;
            this.Density = density;
            this.HeatCapacity = heatCapacity;
        }

        public Material(BsonDocument materialInfo)
        {
            if (materialInfo == null)
            {
                throw new ArgumentNullException(nameof(materialInfo));
            }

            InitializeName(materialInfo.GetValue(MaterialConsts.Name));
            InitializeColor(materialInfo.GetValue(MaterialConsts.Color));
            InitializeThermalConductivity(materialInfo.GetValue(MaterialConsts.ThermalConductivity));
            InitializeDensity(materialInfo.GetValue(MaterialConsts.Density));
            InitializeHeatCapacity(materialInfo.GetValue(MaterialConsts.HeatCapacity));
        }

        public string Name { get; private set; }

        public Color Color { get; private set; }

        public double ThermalConductivity { get; private set; }

        public double Density { get; private set; }

        public double HeatCapacity { get; private set; }

        public BsonDocument ToBsonDocument()
        {
            var document = new BsonDocument
            {
                { MaterialConsts.Name, this.Name },
                { MaterialConsts.Color, $"#{this.Color.R.ToString("X2")}{this.Color.G.ToString("X2")}{this.Color.B.ToString("X2")}" },
                { MaterialConsts.ThermalConductivity, this.ThermalConductivity.ToString() },
                { MaterialConsts.Density, this.Density.ToString() },
                { MaterialConsts.HeatCapacity, this.HeatCapacity.ToString() }
            };

            return document;
        }

        private void InitializeName(BsonValue value)
        {
            Debug.Assert(value != null);

            this.Name = value.ToString();
        }

        private void InitializeColor(BsonValue value)
        {
            Debug.Assert(value != null);

            this.Color = (Color)ColorConverter.ConvertFromString(value.ToString());
        }

        private void InitializeThermalConductivity(BsonValue value)
        {
            Debug.Assert(value != null);

            var correctValue = CorrectDoubleNumber(value.ToString());
            this.ThermalConductivity = double.Parse(correctValue);
        }

        private void InitializeDensity(BsonValue value)
        {
            Debug.Assert(value != null);

            var correctValue = CorrectDoubleNumber(value.ToString());
            this.Density = double.Parse(correctValue);
        }

        private void InitializeHeatCapacity(BsonValue value)
        {
            Debug.Assert(value != null);

            var correctValue = CorrectDoubleNumber(value.ToString());
            this.HeatCapacity = double.Parse(correctValue);
        }

        private string CorrectDoubleNumber(string doubleNumber)
        {
            return doubleNumber.Replace(".", ",");
        }
    }
}
