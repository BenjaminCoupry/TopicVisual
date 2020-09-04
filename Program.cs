using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Drawing;
namespace TopicVisual
{
    static class Program
    {
        
        /// <summary>
        /// Point d'entrée principal de l'application.
        /// </summary>
        [STAThread]
        static void Main()
        {
            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);
            //Application.Run(new Form1());
            Tuple<Base, List<KeyValuePair<string, int>>, List<Tuple<KeyValuePair<string, int>, double[]>>[]> TP = TrouverTopics(@"D:/lab/TopicFinder/TestClusters.txt", 10, 300, 60);
            affTopics(TP,6,0,1,10000,10000,6,100).Save(@"D:/lab/TopicFinder/TestAffClusters1.bmp");
            affTopics(TP, 6, 2, 3, 10000, 10000, 6, 100).Save(@"D:/lab/TopicFinder/TestAffClusters2.bmp");
            affTopics(TP, 6, 4, 5, 10000, 10000, 6, 100).Save(@"D:/lab/TopicFinder/TestAffClusters3.bmp");
        }
        static List<string[]> ExtrairePhrases(string Path)
        {
            List<string[]> retour = new List<string[]>();
            string[] lignes = System.IO.File.ReadAllLines(@Path);
            foreach (string l in lignes)
            {
                string[] phrases = l.Split(new char[] { '.', '!', ';', '(', ')' });
                foreach (string p in phrases)
                {
                    retour.Add(p.Split());
                }
            }
            return retour;
        }
        static List<KeyValuePair<string, int>> ExtraireMots(List<string[]> Phrases)
        {
            Console.WriteLine("Lecture");
            Dictionary<string, int> dict = new Dictionary<string, int>();
            foreach (string[] phrase in Phrases)
            {
                foreach (string mot in phrase)
                {
                    if (!dict.ContainsKey(mot))
                    {
                        dict.Add(mot, 1);
                    }
                    else
                    {
                        dict[mot] += 1;
                    }
                }
            }
            Console.WriteLine("Tri");
            List<KeyValuePair<string, int>> listdict = dict.ToList();
            listdict.Sort((pair1, pair2) => -pair1.Value.CompareTo(pair2.Value));
            return listdict;
        }
        static List<int> PhrasesContenant(List<string[]> Phrases, string Mot)
        {
            List<int> retour = new List<int>();
            for (int i = 0; i < Phrases.Count; i++)
            {
                if (Phrases.ElementAt(i).Contains(Mot))
                {
                    retour.Add(i);
                }
            }
            return retour;
        }
        static int CompterIntersections(List<int> A, List<int> B)
        {
            List<int> courte;
            List<int> longue;
            if (A.Count < B.Count)
            {
                courte = A;
                longue = B;
            }
            else
            {
                courte = B;
                longue = A;
            }
            int compte = 0;
            foreach (int i in courte)
            {
                if (longue.Contains(i))
                {
                    compte++;
                }
            }
            return compte;
        }
        //Renvoie une liste des mots et leur nombre d'occurences, ainsi que la liste des relations entre les mots et les dimsInit mots les plus présents
        static Tuple<List<KeyValuePair<string, int>>, List<double[]>> RelationsMots(List<string[]> Phrases, int dimsInit)
        {
            
            List<KeyValuePair<string, int>> mots = ExtraireMots(Phrases);
            int nbMots = mots.Count;
            //mapping[i] = liste des indices des phrases contenant le mot i
            List<int>[] mapping = new List<int>[nbMots];
            Console.WriteLine("Mapping");
            for (int i = 0; i < nbMots; i++)
            {
                mapping[i] = PhrasesContenant(Phrases, mots.ElementAt(i).Key);
            }
            List<double[]> Relations = new List<double[]>();
            int nombrePhrases = Phrases.Count;
            Console.WriteLine("Construction");
            for (int i = 0; i < nbMots; i++)
            {
                if (i % 100 == 0)
                {
                    Console.WriteLine(i + "/" + nbMots);
                }
                double[] rel = new double[dimsInit];
                for (int j = 0; j < dimsInit; j++)
                {
                    //Nombre de phrases contenant le mot i et le mot j
                    int nbCommun = CompterIntersections(mapping[i], mapping[j]);
                    double pcommun = (double)nbCommun / (double)nombrePhrases;
                    double pi = (double)mapping[i].Count/ (double)nombrePhrases;
                    double pj = (double)mapping[j].Count / (double)nombrePhrases;
                    if (nbCommun != 0)
                    {
                        rel[j] = pcommun/(pi*pj);
                    }
                    else
                    {
                        rel[j] = 0;
                    }
                }
                Relations.Add(rel);
            }
            return new Tuple<List<KeyValuePair<string, int>>, List<double[]>>(mots, Relations);
        }
        static Tuple<List<double[]>, Base> ACPMots(Tuple<List<KeyValuePair<string, int>>, List<double[]>> relMots, int nbComp)
        {
            Random r = new Random();
            return PasserDansBaseAcp(relMots.Item2, 0, nbComp, 0.001, ref r);
        }


        //FONCTION PRINCIPALE
        //Un topic est : Une base dans l'espace des mots, La liste des mots de la base et leurs occurences, un tableau de topic, un topic etant : une liste de mots+nbapparitions avec
        // leur coordonnée dans la base fournie
        static Tuple<Base, List<KeyValuePair<string, int>>, List<Tuple<KeyValuePair<string, int>, double[]>>[]> TrouverTopics(string Path, int nbAcp, int nbMotsDim, int nbTopics)
        {
            Console.WriteLine("Extraction des phrases");
            List<string[]> ph = ExtrairePhrases(Path);
            Console.WriteLine("Calcul des relations");
            Tuple<List<KeyValuePair<string, int>>, List<double[]>> rel = RelationsMots(ph, nbMotsDim);
            Console.WriteLine("Reduction des dim du problème");
            Tuple<List<double[]>, Base> nouvEspaceMots = ACPMots(rel, nbAcp);
            Console.WriteLine("Formation des clusters");
            List<int>[] kmoy = Kmoyennes(nouvEspaceMots.Item1, nbTopics);
            Console.WriteLine("Mise en forme");
            List<Tuple<KeyValuePair<string, int>, double[]>>[] retour = new List<Tuple<KeyValuePair<string, int>, double[]>>[nbTopics];
            for (int i = 0; i < nbTopics; i++)
            {
                List<Tuple<KeyValuePair<string, int>, double[]>> motsAssocies = new List<Tuple<KeyValuePair<string, int>, double[]>>();
                for (int j = 0; j < kmoy[i].Count; j++)
                {
                    int indMot = kmoy[i].ElementAt(j);
                    motsAssocies.Add(new Tuple<KeyValuePair<string, int>, double[]>(rel.Item1.ElementAt(indMot), nouvEspaceMots.Item1.ElementAt(indMot)));
                }
                retour[i] = motsAssocies;
            }
            return new Tuple<Base, List<KeyValuePair<string, int>>, List<Tuple<KeyValuePair<string, int>, double[]>>[]>(nouvEspaceMots.Item2,rel.Item1,retour);
        }
        static void affConsTopics(List<string>[] Topics)
        {
            for (int i = 0; i < Topics.GetLength(0); i++)
            {
                Console.WriteLine("Topic " + i);
                foreach (string st in Topics[i])
                {
                    Console.WriteLine(st);
                }
                Console.WriteLine("------------------------------------------------");
            }
        }
        static Bitmap affTopics(Tuple<Base, List<KeyValuePair<string, int>>, List<Tuple<KeyValuePair<string, int>, double[]>>[]> Topics,int MotsParVect,int compX, int compY, int X, int Y, int policeMin, int policeMax)
        {
            Console.WriteLine("Creation du bitmap");
            Bitmap retour = new Bitmap(X+1, Y+1);
            int nbTop = Topics.Item3.GetLength(0);
            Color[] couleurs = GetNCouleurs(nbTop);
            List<double[]> ech = TopicToEchant(Topics.Item3);
            double[] min = minEchant(ech);
            double[] max = maxEchant(ech);
            int maxApparitionMot = Topics.Item3.Max((pair1) => {
                if (pair1.Count > 0)
                {
                    return pair1.Max((pair2) => pair2.Item1.Value);
                }
                else
                {
                    Console.WriteLine("Classe vide");
                    return 0;
                }
            }
            );
            Tuple<string, List<string>> baseEnMots = BaseEnMots(Topics.Item1, Topics.Item2, MotsParVect);
            int policeN = (policeMax + policeMin) / 2;
            Font drawFont = new Font("Arial", policeN);
            StringFormat drawFormat = new StringFormat();
            drawFormat.FormatFlags = StringFormatFlags.DirectionRightToLeft;
            SolidBrush drawBrush = new SolidBrush(Color.Black);
            Console.WriteLine("Ecriture");
            using (Graphics grf = Graphics.FromImage(retour))
            {
                grf.FillRectangle(drawBrush, 0, 0, X, Y);
                drawBrush = new SolidBrush(Color.White);
                grf.DrawString("Origine : "+baseEnMots.Item1, drawFont, drawBrush, X, policeN, drawFormat);
                string str1 = "Axe y : " + baseEnMots.Item2.ElementAt(compY);
                string str2 = "Axe x : " + baseEnMots.Item2.ElementAt(compX);
                grf.DrawString(str1, drawFont, drawBrush, X, Y-5*policeN, drawFormat);
                grf.DrawString(str2, drawFont, drawBrush, X, Y-2*policeN, drawFormat);
                grf.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.AntiAlias;
                grf.InterpolationMode = System.Drawing.Drawing2D.InterpolationMode.HighQualityBicubic;
                grf.PixelOffsetMode = System.Drawing.Drawing2D.PixelOffsetMode.HighQuality;
                double sx0 = (0 - min[compX]) / (max[compX] - min[compX]);
                double sy0 = (0 - min[compY]) / (max[compY] - min[compY]);
                int x0 = Convert.ToInt32(X * sx0);
                int y0 = Convert.ToInt32(Y * sy0);
                int Rorigine = 5;
                
                for (int i = 0; i < nbTop; i++)
                {
                    Console.WriteLine(i);
                    Color color = couleurs[i];
                    SolidBrush writeBrush = new SolidBrush(color);
                    List<Tuple<KeyValuePair<string, int>, double[]>> Top = Topics.Item3[i];
                    foreach(Tuple<KeyValuePair<string, int>, double[]> mot in Top)
                    {
                        double[] vect = mot.Item2;
                        double sx = (vect[compX] - min[compX]) / (max[compX] - min[compX]);
                        double sy = (vect[compY] - min[compY]) / (max[compY] - min[compY]);
                        int x = Convert.ToInt32(X * sx);
                        int y = Convert.ToInt32(Y * sy);
                        KeyValuePair<string, int> m = mot.Item1;
                        double taille = Math.Log(1.0+ (double)m.Value / (double)maxApparitionMot)/Math.Log(2);
                        int police = (int)((double)policeMin + taille * (double)(policeMax - policeMin));
                        drawFont = new Font("Arial", police);
                        grf.DrawString(m.Key, drawFont, writeBrush, x, y, drawFormat);
                    }
                }
                grf.DrawLine(new Pen(Brushes.White, 1), new Point(x0, y0 - Rorigine), new Point(x0, y0 + Rorigine));
                grf.DrawLine(new Pen(Brushes.White, 1), new Point(x0 - Rorigine, y0), new Point(x0 + Rorigine, y0));
            }
            return retour;
        }


        //Outils couleurs
        static Color[] GetNCouleurs(int N)
        {
            Random r = new Random();
            Color[] retour = new Color[N];
            for(int i=0;i<N;i++)
            {
                retour[i] = Color.FromArgb(r.Next(0, 255), r.Next(0, 255), r.Next(0, 255));
            }
            return retour;
        }
        //Outils diverses
        public static double[] minEchant(List<double[]> echant)
        {
            double[] retour = copiedb(echant.ElementAt(0));
            for (int i = 0; i < echant.Count; i++)
            {
                double[] compare = echant.ElementAt(i);
                for (int k = 0; k < retour.Count(); k++)
                {
                    if (compare[k] < retour[k])
                    {
                        retour[k] = compare[k];
                    }
                }
            }
            return retour;
        }
        public static double[] maxEchant(List<double[]> echant)
        {
            double[] retour = copiedb(echant.ElementAt(0));
            for (int i = 0; i < echant.Count; i++)
            {
                double[] compare = echant.ElementAt(i);
                for (int k = 0; k < retour.Count(); k++)
                {
                    if (compare[k] > retour[k])
                    {
                        retour[k] = compare[k];
                    }
                }
            }
            return retour;
        }
        public static List<double[]> TopicToEchant(List<Tuple<KeyValuePair<string, int>, double[]>>[] Topics)
        {
            List<double[]> retour = new List<double[]>();
            foreach(List<Tuple<KeyValuePair<string, int>, double[]>> top in Topics)
            {
                foreach(Tuple<KeyValuePair<string, int>, double[]> mot in top)
                {
                    retour.Add(mot.Item2);
                }
            }
            return retour;
        }
        //KMoyennes
        static int IndicePointLePlusProche(double[] objet, List<double[]> Points)
        {
            int retour = 0;
            double dist = double.PositiveInfinity;
            for (int i = 0; i < Points.Count; i++)
            {
                double d = normeVect(SommmeVect(objet, Points.ElementAt(i), true));
                if (d < dist)
                {
                    dist = d;
                    retour = i;
                }
            }
            return retour;
        }
        static double[] PosMoyenne(List<double[]> Objets, List<int> Indices)
        {
            double[] posMoy = new double[Objets.ElementAt(0).Count()];
            for (int i = 0; i < Objets.ElementAt(0).Count(); i++)
            {
                posMoy[i] = 0;
            }
            for (int i = 0; i < Indices.Count; i++)
            {
                posMoy = SommmeVect(posMoy, Objets.ElementAt(Indices.ElementAt(i)), false);
            }
            posMoy = multiplierVect(posMoy, 1.0 / Indices.Count);
            return posMoy;
        }
        static List<double[]> SamplerNPositions(List<double[]> Objets, int n)
        {
            List<double[]> retour = new List<double[]>();
            Random r = new Random();
            List<int> IndicesPossibles = new List<int>();
            for (int i = 0; i < Objets.Count; i++)
            {
                IndicesPossibles.Add(i);
            }
            for (int i = 0; i < n; i++)
            {
                int u = r.Next(0, IndicesPossibles.Count);
                int ind = IndicesPossibles.ElementAt(u);
                IndicesPossibles.RemoveAt(u);
                retour.Add(copiedb(Objets.ElementAt(ind)));
            }
            return retour;
        }
        static bool ClusterDiff(List<double[]> clA, List<double[]> clB, double eps)
        {
            for (int i = 0; i < clB.Count; i++)
            {
                if (normeVect(SommmeVect(clA.ElementAt(i), clB.ElementAt(i), true)) > eps)
                {
                    return true;
                }
            }
            return false;
        }
        static List<int>[] Kmoyennes(List<double[]> Objets, int Kclusters)
        {
            List<double[]> clustersCoord = SamplerNPositions(Objets, Kclusters);
            List<int>[] clusters = new List<int>[Kclusters];
            bool change = true;
            while (change)
            {
                for (int k = 0; k < Kclusters; k++)
                {
                    clusters[k] = new List<int>();
                }
                for (int i = 0; i < Objets.Count; i++)
                {
                    clusters[IndicePointLePlusProche(Objets.ElementAt(i), clustersCoord)].Add(i);
                }
                List<double[]> clustersCoordN = new List<double[]>();
                for (int k = 0; k < Kclusters; k++)
                {
                    clustersCoordN.Add(PosMoyenne(Objets, clusters[k]));
                }
                change = ClusterDiff(clustersCoordN, clustersCoord, 1e-6);
                clustersCoord = clustersCoordN;
            }
            return clusters;
        }

        //base et mots
        static Tuple<string,List<string>> BaseEnMots(Base B, List<KeyValuePair<string, int>> Mots, int motsParVect)
        {
            Console.WriteLine("Transcription de la base ACP en mots");
            List<string> retourL = new List<string>();
            for(int i=0;i<B.vecteurs.Count();i++)
            {
                retourL.Add(VectBaseEnMot(B.vecteurs.ElementAt(i), Mots, motsParVect));
            }
            return new Tuple<string, List<string>>(VectBaseEnMot(B.origine, Mots, motsParVect),retourL);
        }
        static string VectBaseEnMot(double[] vect, List<KeyValuePair<string, int>> Mots, int nbCompAffichees)
        {
            List<Tuple<double, int>> vectTriable = new List<Tuple<double, int>>();
            for(int i=0;i< vect.GetLength(0); i++)
            {
                vectTriable.Add(new Tuple<double, int>(vect[i], i));
            }
            vectTriable.Sort((pair1, pair2) => -Math.Abs(pair1.Item1).CompareTo(Math.Abs(pair2.Item1)));
            double ccordMax = vect.Max((I) => Math.Abs(I));
            string retour = "[";
            for(int i=0; i<nbCompAffichees;i++)
            {
                Tuple<double, int> comp = vectTriable.ElementAt(i);
                KeyValuePair<string, int> mot = Mots.ElementAt(comp.Item2);
                double prevalence = comp.Item1 / ccordMax;
                retour = retour + Math.Abs(prevalence*100).ToString("0.00")+"*"+ mot.Key+"("+mot.Value+")"+";";
            }
            retour = retour + "...]";
            return retour;
        }

        //ACP
        public static Tuple<List<double[]>,Base> PasserDansBaseAcp(List<double[]> echants, int nbParametresStables, int nbCompAcp, double epsilonAcp, ref Random r)
        {
            //Passe les parametres dans la base ACP generee par le set d'echantillons, en gardant la correspondance avec les elements non parametres des echantillons
            List<double[]> echants_tronque = new List<double[]>();
            List<double[]> echants_retour = new List<double[]>();
            Console.WriteLine("Conditionement des données");
            for (int i = 0; i < echants.Count; i++)
            {
                echants_tronque.Add(SousVecteur(echants.ElementAt(i), nbParametresStables, echants.ElementAt(i).Count() - nbParametresStables));
            }
            Console.WriteLine("Calcul de la base");
            Base nouvelleBase = BaseAcp(echants_tronque, nbCompAcp, epsilonAcp, ref r);
            Console.WriteLine("Reconditionement des données");
            for (int i = 0; i < echants.Count; i++)
            {
                double[] vectRetour = SousVecteur(echants.ElementAt(i), 0, nbParametresStables).Concat(CoordoneesDansBaseAcp(nouvelleBase, echants_tronque.ElementAt(i))).ToArray();
                echants_retour.Add(vectRetour);
            }
            return new Tuple<List<double[]>, Base> (echants_retour,nouvelleBase);
        }
        public static Base BaseAcp(List<double[]> echants, int n, double epsilon, ref Random r)
        {
            //Retourne la liste des vecteurs de la base et le vecteur moyen, par ordre decroissant
            List<double> valeursPropres = new List<double>();
            Console.WriteLine("Centrer les donées...");
            Tuple<double[,], double[]> EchantEtCentre = MatriceEchantsCentres(echants);
            double[] moy = EchantEtCentre.Item2;
            List<Tuple<double[], double>> vp;
            List<double[]> vect_base = new List<double[]>();
            Console.WriteLine("Calcul des vp...");
            vp = ACP(EchantEtCentre.Item1, n, epsilon, ref r);
            Console.WriteLine("Export...");
            for (int i = 0; i < n; i++)
            {
                double[] v = vp.ElementAt(i).Item1;
                vect_base.Add(v);
                valeursPropres.Add(vp.ElementAt(i).Item2);
                Console.WriteLine(i + " : " + vp.ElementAt(i).Item2);
            }

            return new Base(vect_base, valeursPropres, moy);
        }
        public static double[] CoordoneesDansBaseAcp(Base Baseacp, double[] vecteur)
        {
            double[] V = SommmeVect(vecteur, Baseacp.origine, true);
            List<double> retour = new List<double>();
            for (int i = 0; i < Baseacp.vecteurs.Count; i++)
            {
                retour.Add(produitScalaire(V, Baseacp.vecteurs.ElementAt(i)));
            }
            return retour.ToArray();
        }
        public static double[] SousVecteur(double[] V, int n, int m)
        {
            //Retourne un sous vecteur pribvé de ses n premieres coordonnées
            return V.ToList().GetRange(n, m).ToArray();
        }
        public static List<Tuple<double[], double>> ACP(double[,] Mc, int n, double epsilon, ref Random r)
        {
            //Acp pour petit nombre de variables a decorreler
            double K = Mc.GetLength(0);
            Console.WriteLine("Calcul des var covar");
            double[,] varcovar = ScalaireMat(1.0 / K, ProduitMat(Transpose(Mc), Mc));
            List<Tuple<double[], double>> retour = new List<Tuple<double[], double>>();
            for (int i = 0; i < n; i++)
            {
                Console.WriteLine("vp " + i + "/" + n);
                Tuple<double[], double> vp = PuissanceIteree(varcovar, VectAlea(varcovar.GetLength(0), ref r), epsilon);
                varcovar = Deflation(varcovar, vp);
                retour.Add(vp);
            }
            return retour;
        }
        private static double[,] EchantsVersMat(ref List<double[]> echants)
        {
            //Renvoie une matrice dont chaque ligne est un echantillon
            int K = echants.Count();
            int N = echants.ElementAt(0).GetLength(0);
            double[,] M = new double[K, N];
            for (int k = 0; k < K; k++)
            {
                double[] v = echants.ElementAt(k);
                for (int n = 0; n < N; n++)
                {
                    M[k, n] = v[n];
                }
            }
            return M;
        }
        private static double[,] Centrer(double[,] M)
        {
            Console.WriteLine(">Centrer");
            int l = M.GetLength(0);
            int c = M.GetLength(1);
            double[,] retour = new double[l, c];
            for (int j = 0; j < c; j++)
            {
                if (j % 100 == 0)
                {
                    Console.WriteLine(j + "/" + c);
                }
                double[] var = MatVersVect(M, j);
                double Esp = var.Average();
                for (int i = 0; i < l; i++)
                {
                    retour[i, j] = M[i, j] - Esp;
                }
            }
            return retour;
        }
        private static double[] IndivMoyen(double[,] M)
        {
            Console.WriteLine(">moyenne");
            //Renvoie l'individu moyen a partir de la matrice des echantillons
            int n = M.GetLength(1);
            double[] retour = new double[n];
            for (int i = 0; i < n; i++)
            {
                if (i % 100 == 0)
                {
                    Console.WriteLine(i + "/" + n);
                }
                retour[i] = MatVersVect(M, i).Average();
            }
            return retour;
        }
        private static Tuple<double[,], double[]> MatriceEchantsCentres(List<double[]> echants)
        {
            //Renvoie la matrice centree des echantillons et sa moyenne
            double[,] M = EchantsVersMat(ref echants);
            double[] C = IndivMoyen(M);
            return new Tuple<double[,], double[]>(Centrer(M), C);
        }
        public static Tuple<double[], double> PuissanceIteree(double[,] M, double[] v0, double epsilon)
        {
            double[] x = copiedb(v0);
            double beta = ProduitMat(Transpose(x), ProduitMat(M, x))[0, 0];
            double beta_t;
            do
            {
                x = normerVect(MatVersVect(ProduitMat(M, x), 0));
                beta_t = beta;
                beta = ProduitMat(Transpose(x), ProduitMat(M, x))[0, 0];
            }
            while (Math.Abs((beta - beta_t) / beta_t) > epsilon);
            double s = 1.0;
            if (x[0] < 0)
            {
                s = -1.0;
            }
            return new Tuple<double[], double>(multiplierVect(x, s), beta);
        }
        public static double[,] Deflation(double[,] M, Tuple<double[], double> vp)
        {
            return SommeMat(M, ScalaireMat(-vp.Item2, ProduitMat(vp.Item1, Transpose(vp.Item1))));
        }
        private static double[] VectAlea(int n, ref Random r)
        {
            double[] vect = new double[n];
            for (int i = 0; i < n; i++)
            {
                vect[i] = 2.0 * r.NextDouble() - 1.0;
            }
            return vect;
        }
        public class Base
        {
            public List<double[]> vecteurs;
            public List<double> valeursPropres;
            public double[] origine;

            public Base(List<double[]> vecteurs, List<double> valeursPropres, double[] origine)
            {
                this.vecteurs = vecteurs;
                this.valeursPropres = valeursPropres;
                this.origine = origine;
            }
        }

        private static double produitScalaire(double[] V1, double[] V2)
        {
            //Produit scalaire de deux vecteurs
            int Taille = Math.Min(V1.Length, V2.Length);
            double Somme = 0.0;
            for (int i = 0; i < Taille; i++)
            {
                Somme += V1[i] * V2[i];
            }
            return Somme;
        }
        public static double[] SommmeVect(double[] a, double[] b, bool sub)
        {
            double[] retour = new double[a.GetLength(0)];
            double sgn = 1.0;
            if (sub)
            {
                sgn = -1.0;
            }
            for (int i = 0; i < a.GetLength(0); i++)
            {
                retour[i] = a[i] + sgn * b[i];
            }
            return retour;
        }
        public static double[,] Transpose(double[,] M)
        {
            //Matrice transposee
            double[,] retour = new double[M.GetLength(1), M.GetLength(0)];
            for (int i = 0; i < M.GetLength(1); i++)
            {
                for (int j = 0; j < M.GetLength(0); j++)
                {
                    retour[i, j] = M[j, i];
                }
            }
            return retour;
        }
        public static double[,] Transpose(double[] M)
        {
            return Transpose(VectVersMat(M));
        }
        public static double[,] ScalaireMat(double k, double[,] M)
        {
            //Multiplie la matrice par un scalaire
            double[,] retour = new double[M.GetLength(0), M.GetLength(1)];
            for (int i = 0; i < M.GetLength(0); i++)
            {
                for (int j = 0; j < M.GetLength(1); j++)
                {
                    retour[i, j] = k * M[i, j];
                }
            }
            return retour;
        }
        public static double[,] ProduitMat(double[,] A, double[,] B)
        {
            //Produit matriciel
            double[,] retour = new double[A.GetLength(0), B.GetLength(1)];
            for (int i = 0; i < A.GetLength(0); i++)
            {
                for (int j = 0; j < B.GetLength(1); j++)
                {
                    double somme = 0;
                    for (int k = 0; k < A.GetLength(1); k++)
                    {
                        somme += A[i, k] * B[k, j];
                    }
                    retour[i, j] = somme;
                }
            }
            return retour;
        }
        public static double[,] ProduitMat(double[,] A, double[] B)
        {
            return ProduitMat(A, VectVersMat(B));
        }
        public static double[,] ProduitMat(double[] A, double[,] B)
        {
            return ProduitMat(VectVersMat(A), B);
        }
        public static double[,] VectVersMat(double[] v)
        {
            //transforme un vecteur en matrice
            double[,] retour = new double[v.Length, 1];
            for (int i = 0; i < v.Length; i++)
            {
                retour[i, 0] = v[i];
            }
            return retour;
        }
        public static double[] MatVersVect(double[,] M, int p)
        {
            //transforme une matrice en vecteur
            double[] retour = new double[M.GetLength(0)];
            for (int i = 0; i < M.GetLength(0); i++)
            {
                retour[i] = M[i, p];
            }
            return retour;
        }
        public static double[,] SommeMat(double[,] A, double[,] B)
        {
            //Somme des matrices
            double[,] retour = new double[A.GetLength(0), A.GetLength(1)];
            for (int i = 0; i < A.GetLength(0); i++)
            {
                for (int j = 0; j < A.GetLength(1); j++)
                {
                    retour[i, j] = A[i, j] + B[i, j];
                }
            }
            return retour;
        }
        private static double[] multiplierVect(double[] V, double k)
        {
            //Multiplier le vecteur par k
            int Taille = V.Length;
            double[] ret = new double[Taille];
            for (int i = 0; i < Taille; i++)
            {
                ret[i] = V[i] * k;
            }
            return ret;
        }
        public static double[] normerVect(double[] u)
        {
            double[] un = copiedb(u);
            multiplierVect(ref un, 1.0 / normeVect(un));
            return un;
        }
        private static double normeVect(double[] V)
        {
            //Norme du vecteur v
            return Math.Sqrt(produitScalaire(V, V));
        }
        private static void multiplierVect(ref double[] V, double k)
        {
            //Multiplier le vecteur par k
            int Taille = V.Length;
            for (int i = 0; i < Taille; i++)
            {
                V[i] *= k;
            }
        }
        private static double[] copiedb(double[] V1)
        {
            //renvoie une copie de V1
            int Taille = V1.Length;
            double[] Resultat = new double[Taille];
            for (int i = 0; i < Taille; i++)
            {
                Resultat[i] = V1[i];
            }
            return Resultat;
        }

    }
}
