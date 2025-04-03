var version,header;
function autoversion()
{
 this.active=false;
};
autoversion.prototype.includeFile=function(fn)
{
 var fa,ver="";
 if(this.active && header && header.version){ver="_ver"+header.version;};
 if(fn.indexOf('.js')!=-1)
 {
  fa=fn.split('.js');
  document.write("<script type='text/javascript' src='"+fa[0]+ver+".js"+fa[1]+"'><\/script>\n");
  return;
 };
 if(fn.indexOf('.css')!=-1)
 {
  fa=fn.split('.css');
  document.write("<link rel='stylesheet' type='text/css' href='"+fa[0]+ver+".css"+fa[1]+"'>\n");
  return;
 };
 alert('cannot include '+fn);
};
autoversion.prototype.includeHeader=function()
{
 if(this.active)
 {
  this.includeFile("header_ver"+(Math.floor(Math.random() * 1000) + 1)+".js");
 }
 else
 {
  this.includeFile("header.js?ver="+(Math.floor(Math.random() * 1000) + 1));
 };
};
autoversion.prototype.includeFiles=function(fa,v)
{
 var i;
 for(i=0;i<fa.length;i++){this.include_file(fn);};
};
version=new autoversion();
version.includeFile('../fs/fs_ver_ver000.js');
