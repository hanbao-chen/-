"""Build the signed PDF from the reviewed Markdown; local images stay local."""
from pathlib import Path
import argparse,re,html,json
from reportlab.pdfgen import canvas
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,Image,PageBreak,KeepTogether
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT,TA_CENTER,TA_JUSTIFY
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

ROOT=Path(__file__).resolve().parents[1]
REPORT=ROOT/'report'

def main(body_font,bold_font):
    pdfmetrics.registerFont(TTFont('BodyCJK',str(body_font),subfontIndex=0))
    pdfmetrics.registerFont(TTFont('BoldCJK',str(bold_font),subfontIndex=0))
    W,H=A4;usable=W-96
    styles={
        'title':ParagraphStyle('Title',fontName='BoldCJK',fontSize=21,leading=30,spaceAfter=15,wordWrap='CJK'),
        'h2':ParagraphStyle('Heading',fontName='BoldCJK',fontSize=13.4,leading=20,spaceBefore=4,spaceAfter=10,wordWrap='CJK',keepWithNext=True),
        'h3':ParagraphStyle('Subheading',fontName='BoldCJK',fontSize=11,leading=17,spaceBefore=5,spaceAfter=6,wordWrap='CJK',keepWithNext=True),
        'body':ParagraphStyle('Body',fontName='BodyCJK',fontSize=10.2,leading=16.3,spaceAfter=8,wordWrap='CJK',splitLongWords=True),
        'caption':ParagraphStyle('Caption',fontName='BodyCJK',fontSize=8.3,leading=12.5,spaceAfter=8,wordWrap='CJK'),
        'table':ParagraphStyle('Table',fontName='BodyCJK',fontSize=8.5,leading=12.6,wordWrap='CJK'),
        'table_head':ParagraphStyle('TableHead',fontName='BoldCJK',fontSize=8.5,leading=12.6,wordWrap='CJK'),
        'code':ParagraphStyle('Code',fontName='BodyCJK',fontSize=8,leading=12.5,spaceAfter=3,wordWrap='CJK'),
    }
    def markup(s):
        s=html.escape(s)
        s=re.sub(r'\[([^\]]+)\]\(([^)]+)\)',lambda m:f'<link href="{html.escape(html.unescape(m.group(2)),quote=True)}" color="#15537c">{m.group(1)}</link>',s)
        s=re.sub(r'\*\*(.+?)\*\*',r'<font name="BoldCJK">\1</font>',s)
        s=re.sub(r'`([^`]+)`',r'\1',s)
        return s
    def para(s,style='body'):return Paragraph(markup(s),styles[style])
    source=(REPORT/'course_report_22634107.md').read_text('utf-8')
    lines=source.splitlines();story=[];i=0;in_code=False
    while i<len(lines):
        line=lines[i].strip()
        if line=='<!-- PAGE -->':story.append(PageBreak());i+=1;continue
        if line.startswith('```'):in_code=not in_code;i+=1;continue
        if in_code:
            if line:story.append(para(line,'code'))
            i+=1;continue
        if not line:i+=1;continue
        if line.startswith('# '):story.append(para(line[2:],'title'));i+=1;continue
        if line.startswith('### '):story.append(para(line[4:],'h3'));i+=1;continue
        if line.startswith('## '):story.append(para(line[3:],'h2'));i+=1;continue
        if line.startswith('|'):
            rows=[]
            while i<len(lines) and lines[i].strip().startswith('|'):
                cells=[s.strip() for s in lines[i].strip().strip('|').split('|')]
                if not all(re.fullmatch(r':?-+:?',s or '-') for s in cells):rows.append(cells)
                i+=1
            cols=len(rows[0]);widths=([usable*.32,usable*.68] if cols==2 else ([usable*.46]+[usable*.54/(cols-1)]*(cols-1)))
            data=[[para(s,'table_head' if r==0 else 'table') for s in row] for r,row in enumerate(rows)]
            table=Table(data,colWidths=widths,hAlign='LEFT',repeatRows=1)
            table.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e7edf2')),('GRID',(0,0),(-1,-1),.4,colors.HexColor('#cbd2d8')),('VALIGN',(0,0),(-1,-1),'MIDDLE'),('LEFTPADDING',(0,0),(-1,-1),7),('RIGHTPADDING',(0,0),(-1,-1),7),('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f7f9fb')])]))
            story.extend([table,Spacer(1,10)]);continue
        m=re.match(r'!\[([^\]]*)\]\(([^)]+)\)',line)
        if m:
            p=(REPORT/m.group(2)).resolve();im=Image(str(p));scale=min(usable/im.imageWidth,285/im.imageHeight)
            im.drawWidth=im.imageWidth*scale;im.drawHeight=im.imageHeight*scale;im.hAlign='CENTER'
            group=[im,Spacer(1,5)]
            j=i+1
            while j<len(lines) and not lines[j].strip():j+=1
            if j<len(lines) and lines[j].strip().startswith('图'):
                group.append(para(lines[j].strip(),'caption'));i=j
            story.append(KeepTogether(group));i+=1;continue
        style='caption' if line.startswith('图') else 'body'
        story.append(para(line,style));i+=1
    def footer(c,doc):
        c.saveState();c.setFont('BodyCJK',8);c.setFillColor(colors.HexColor('#555555'))
        c.drawString(48,25,'海洋数据挖掘技术  朱宸浩  22634107')
        c.drawRightString(W-48,25,str(doc.page));c.restoreState()
    pdf=REPORT/'course_report_22634107.pdf'
    doc=SimpleDocTemplate(str(pdf),pagesize=A4,rightMargin=48,leftMargin=48,topMargin=43,bottomMargin=43,title='近岸 H3 小海域热点识别中的证据与代价',author='朱宸浩 22634107',subject='海洋数据挖掘技术课程研究报告',allowSplitting=True)
    doc.build(story,onFirstPage=footer,onLaterPages=footer)
    print(pdf)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--body-font',type=Path,default=Path(r'C:\Windows\Fonts\simsun.ttc'));p.add_argument('--bold-font',type=Path,default=Path(r'C:\Windows\Fonts\simhei.ttf'));a=p.parse_args();main(a.body_font,a.bold_font)
