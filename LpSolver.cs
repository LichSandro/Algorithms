using System;
using System.Globalization;
using System.Linq;
using System.Threading;

namespace LinearProgramming
{
	public class LpSolver
	{
		private const double Eps = 1.0e-9;

		// returns max c*x such that A*x <= b, x >= 0
		public static double? SolveStandardLpProblem(double[][] A, double[] b, double[] c, out double[] x)
		{
			int m = A.Length;
			int n = A[0].Length + 1;
			x = new double[n - 1];

			var index = new int[n + m];
			for (int i = 0; i < n + m; i++)
				index[i] = i;

			var a = new double[m + 2][];
			for (var i = 0; i < a.Length; i++)
				a[i] = Enumerable.Repeat(0.0, n + 1).ToArray();

			int L = m;
			for (int i = 0; i < m; i++)
			{
				for (int j = 0; j < n - 1; j++)
					a[i][j] = -A[i][j];

				a[i][n - 1] = 1.0;
				a[i][n] = b[i];
				if (a[L][n] > a[i][n])
					L = i;
			}

			for (int j = 0; j < n - 1; j++)
				a[m][j] = c[j];

			a[m + 1][n - 1] = -1.0;
			for (int E = n - 1;;)
			{
				if (L < m)
				{
					int t = index[E];
					index[E] = index[L + n];
					index[L + n] = t;
					a[L][E] = 1.0 / a[L][E];
					for (int j = 0; j <= n; j++)
					{
						if (j != E)
							a[L][j] = -a[L][j] * a[L][E];
					}

					for (int i = 0; i <= m + 1; i++)
					{
						if (i != L)
						{
							for (int j = 0; j <= n; j++)
							{
								if (j != E)
									a[i][j] += a[L][j] * a[i][E];
							}

							a[i][E] = a[i][E] * a[L][E];
						}
					}
				}

				E = -1;
				for (int j = 0; j < n; j++)
				{
					if (E < 0 || index[E] > index[j])
					{
						if (a[m + 1][j] > Eps || (Math.Abs(a[m + 1][j]) < Eps && a[m][j] > Eps))
							E = j;
					}
				}

				if (E < 0)
					break;

				L = -1;
				for (int i = 0; i < m; i++)
				{
					if (a[i][E] < -Eps)
					{
						if (L < 0)
							L = i;
						else
						{
							var d = a[L][n] / a[L][E] - a[i][n] / a[i][E];
							if (d < -Eps || (Math.Abs(d) < Eps && index[L + n] > index[i + n]))
								L = i;
						}
					}
				}

				if (L < 0)
					return double.PositiveInfinity;
			}

			if (a[m + 1][n] < -Eps)
				return null;

			for (int i = 0; i < m; i++)
				if (index[n + i] < n - 1)
					x[index[n + i]] = a[i][n];

			return a[m][n];
		}
	}

	public class Program
	{
		public static void SolveLpProblem(double[][] a, double[] b, double[] c)
		{
			var res = LpSolver.SolveStandardLpProblem(a, b, c, out var x);
			if (res == null)
			{
				Console.WriteLine("No feasible solution exists, the constraints are inconsistent!");
				return;
			}

			if (double.IsPositiveInfinity(res.Value))
			{
				Console.WriteLine("Answer = +Inf, the polytope is unbounded in the direction of the gradient of the objective function!");
				return;
			}

			Console.WriteLine($"Answer = {res}, x = [{string.Join(", ", x)}]");
		}

		public static void Main()
		{
			Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;

			double[][] a1 =
			{
				new[] {1.0, -1.0},
				new[] {-1.0, 1.0},
				new[] {1.0, 1.0},
			};

			double[] b1 = {1.0, 1.0, 3.0};

			double[] c1 = {2.0, 1.0};

			SolveLpProblem(a1, b1, c1);

			Console.WriteLine();

			double[][] a2 =
			{
				new[] {1.0, -1.0},
				new[] {-1.0, 1.0},
				new[] {-1.0, -1.0},
			};

			double[] b2 = {1.0, 1.0, -3.0};

			double[] c2 = {2.0, 1.0};

			SolveLpProblem(a2, b2, c2);

			Console.WriteLine();

			double[][] a3 =
			{
				new[] {0.0, -1.0},
				new[] {-1.0, 0.0},
				new[] {1.0, 1.0},
			};

			double[] b3 = {-2.0, -2.0, 1.0};

			double[] c3 = {2.0, 1.0};

			SolveLpProblem(a3, b3, c3);
		}
	}
}