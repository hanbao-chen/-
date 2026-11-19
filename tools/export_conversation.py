"""Export only actual user/assistant messages and structured user questions."""
from pathlib import Path
import json,argparse,hashlib,re,shutil
from datetime import datetime,timezone,timedelta

ROOT=Path(__file__).resolve().parents[1]

def export(session):
    items=[]
    for line in session.read_text('utf-8').splitlines():
        try:items.append(json.loads(line))
        except json.JSONDecodeError:continue
    output=['# 海洋数据挖掘技术课程作业完整研究对话','',
        '作者：朱宸浩　学号：22634107','',
        '本文件直接转换本作业会话的实际持久化记录。用户原文、AI对外文字和结构化问题均保留；没有将摘要扩写为虚构对话。',
        '应用自动注入的环境配置、系统与开发者指令、模型内部推理不属于人机研究对话，不导出。代码、实验原始计数与工程核查见仓库experiments目录。',
        '导出时仍在准备和发布交付物；截止时间之后的上传完成通知不在本次快照中。导出时间、记录条数和源记录指纹见export_metadata.json。','']
    nuser=nai=nq=0;last='';user_checks=[];attachments=[]
    for event in items:
        if event.get('type')!='response_item':continue
        p=event.get('payload',{});ts=event.get('timestamp','')
        if p.get('type')=='message' and p.get('role') in ['user','assistant']:
            text='\n'.join(c.get('text','') for c in p.get('content',[]) if c.get('type') in ['input_text','output_text','text'])
            if not text.strip():continue
            if p['role']=='user' and text.lstrip().startswith('<environment_context>'):continue
            role='用户' if p['role']=='user' else 'AI'
            if role=='用户':nuser+=1;user_checks.append(text)
            else:nai+=1
            try:date=datetime.fromisoformat(ts.replace('Z','+00:00')).astimezone(timezone(timedelta(hours=8))).isoformat()
            except ValueError:date=ts
            output.extend([f'## {role} · {date}',''])
            if text.lstrip().startswith('<send_user_message_question_reply>'):
                output.extend(['应用传入的结构化答复原文：','','```xml',text.rstrip(),'```',''])
            else:output.extend([text.rstrip(),''])
            if role=='用户' and any(c.get('type')=='input_image' for c in p.get('content',[])):
                matches=re.findall(r'<image[^>]*path="([^"]+)"',text)
                for original in matches:
                    image=Path(original.replace('\\\\','\\'))
                    if image.is_file():
                        folder=ROOT/'conversations'/'attachments';folder.mkdir(parents=True,exist_ok=True)
                        name=f'attachment_{len(attachments)+1:03d}'+image.suffix.lower()
                        shutil.copy2(image,folder/name)
                        relative='attachments/'+name
                        attachments.append(relative)
                        output.extend(['![用户上传的原始截图]('+relative+')',''])
            last=ts
        elif p.get('type') in ['function_call','custom_tool_call'] and 'request_user_input_async' in p.get('name',''):
            try:args=json.loads(p.get('arguments','{}'))
            except (json.JSONDecodeError,TypeError):args={}
            for q in args.get('questions',[]):
                nq+=1;output.extend([f'## AI结构化提问 · {ts}','',q.get('title','')]+['- '+s for s in q.get('options',[])]+[''])
    dest=ROOT/'conversations';dest.mkdir(exist_ok=True)
    text='\n'.join(output)
    required=['其实没有理解','先理解独立观测','理解看起来更好的错误结果','弄清统计数字的含义','讨论严格校正的代价','把等待联系到实际应用','追问门槛附近的困难','检查“待观察”的实际价值','确定作业的贡献','明确结论边界','形成自己的理解','22634107']
    checked={s:s in text for s in required}
    if not all(checked.values()):raise RuntimeError('Research dialogue completeness check failed: '+repr(checked))
    (dest/'full_conversation.md').write_text(text,encoding='utf-8')
    metadata={'exported_at_hong_kong':datetime.now(timezone(timedelta(hours=8))).isoformat(),'source_session':session.name,'source_bytes':session.stat().st_size,'source_sha256':hashlib.sha256(session.read_bytes()).hexdigest(),'last_exported_message_utc':last,'user_message_count':nuser,'assistant_message_count':nai,'structured_question_count':nq,'image_attachments':attachments,'all_required_research_questions_found':checked,'transcript_sha256':hashlib.sha256((dest/'full_conversation.md').read_bytes()).hexdigest(),'scope':'actual visible research and preparation messages through the export cutoff; delivery notice after export excluded','reconstructed_or_fictional_dialogue':False}
    (dest/'export_metadata.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:metadata[k] for k in ['exported_at_hong_kong','user_message_count','assistant_message_count','structured_question_count','last_exported_message_utc']},ensure_ascii=False,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--session',type=Path,required=True);a=p.parse_args();export(a.session)
