package jp;

import java.io.*;

/****************************/
/* 主成分分析               */
/*      coded by Y.Suganuma */
/****************************/
public class Test {
	public static void main(String args[]) throws IOException
	{
		double x[][], r[], a[][];
		int i1, i2, n, p, sw;
		String str;
		Data dt;
		BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
					// 変数の数とデータの数
		str = in.readLine();
		dt  = new Data(str);
		dt.next();
		p = Integer.parseInt(dt.data);
		dt.next();
		n = Integer.parseInt(dt.data);

		r = new double [p];
		x = new double [p][n];
		a = new double [p][p];
					// データ
		for (i1 = 0; i1 < n; i1++) {
			str = in.readLine();
			dt  = new Data(str);
			for (i2 = 0; i2 < p; i2++) {
				dt.next();
				x[i2][i1] = Double.parseDouble(dt.data);
			}
		}

		sw = App.principal(p, n, x, r, a, 1.0e-10, 200);

		if (sw == 0) {
			for (i1 = 0; i1 < p; i1++) {
				System.out.print("主成分 " + r[i1]);
				System.out.print(" 係数");
				for (i2 = 0; i2 < p; i2++)
					System.out.print(" " + a[i1][i2]);
				System.out.println();
			}
		}
		else
			System.out.println("***error  解を求めることができませんでした");
	}
}