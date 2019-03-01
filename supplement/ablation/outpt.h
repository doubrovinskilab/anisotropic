int q = 10000;
if( !(i/double(q)-i/q) )
{
	char str1[80]="outpt/out";
	char str3[80]=".m";

	strncat(str1,itoa(i/q,10),10);
	strncat(str1,str3,10);

	ofstream * file = new ofstream(str1);

	// file->precision(20);

	if( NumPointsInHole )
	{
		(*file)<<"PointsInHole=[";
		for(j=0; j<NumPointsInHole; j++)
		{
			(*file)<<PointsInHole[j]<<" ";
		}
		(*file)<<"];"<<endl;
	}

	(*file)<<"x=[";
	for(j=0; j<NumPoints; j++)
	{
		(*file)<<x[j]<<" ";
	}
	(*file)<<"];"<<endl;

	(*file)<<"y=[";
	for(j=0; j<NumPoints; j++)
	{
		(*file)<<y[j]<<" ";
	}
	(*file)<<"];"<<endl;

	file->close();
	delete file;
}
