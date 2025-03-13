#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script para enviar e-mails para usuários com informações de acesso ao servidor,
baixando as informações de uma planilha do Google Drive e atualizando o status na planilha.

Autor: Raphael Luiz Lobo da Silva Souza
Data: 13/03/2025
Descrição: Este script autentica no Google Drive e Google Sheets, baixa a planilha com dados de usuários,
envia os e-mails com informações de acesso e atualiza o status de envio na planilha.
"""

import os
import smtplib
import pandas as pd
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from google.oauth2.credentials import Credentials
from google_auth_oauthlib.flow import InstalledAppFlow
from googleapiclient.discovery import build
from googleapiclient.http import MediaIoBaseDownload
import io

# Função para autenticar no Google Drive e Google Sheets
def authenticate_google_drive():
    """Autentica e retorna a API do Google Drive e Google Sheets."""
    creds = None
    SCOPES = ['https://www.googleapis.com/auth/drive.readonly', 'https://www.googleapis.com/auth/spreadsheets']

    if os.path.exists('token.json'):
        creds = Credentials.from_authorized_user_file('token.json', SCOPES)

    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
        else:
            flow = InstalledAppFlow.from_client_secrets_file(
                'credentials.json', SCOPES)
            creds = flow.run_local_server(port=8080)  # A porta fixa é 8080

        with open('token.json', 'w') as token:
            token.write(creds.to_json())

    service = build('drive', 'v3', credentials=creds)
    sheets_service = build('sheets', 'v4', credentials=creds)
    return service, sheets_service

# Função para baixar e exportar o arquivo Excel do Google Drive
def download_file(service, file_id):
    """Baixa e exporta o arquivo do Google Drive (Google Sheets) para .xlsx."""
    request = service.files().export_media(fileId=file_id, mimeType='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
    fh = io.FileIO('usuarios.xlsx', 'wb')  # Salvar como .xlsx
    downloader = MediaIoBaseDownload(fh, request)
    done = False
    while done is False:
        status, done = downloader.next_chunk()
        print(f"Download {int(status.progress() * 100)}%.")
    print("Arquivo Excel baixado com sucesso!")

# Função para enviar e-mail
def send_email(nome, email, senha, usuario, sheets_service):
    from_email = 'raphaluizlobo@gmail.com'
    from_password = 'lzdkdllhbbsxjqqn'  # Substitua pela senha de aplicativo criada

    smtp_server = 'smtp.gmail.com'
    smtp_port = 587

    corpo_email = f"""
    Caro(a) {nome},

    Sua conta no servidor do CEPID B3 foi criada com sucesso. Abaixo seguem as informações necessárias para que você possa acessar o servidor.

    Para acessar, utilize um cliente de SSH e execute o seguinte comando no terminal:

    ssh {usuario}@davinci.icb.usp.br

    Servidor: davinci.icb.usp.br
    Porta: 22
    Usuário: {usuario}
    Senha: {senha}

    ATENÇÃO: A mudança de senha é obrigatória no seu primeiro acesso. Isso pode ser feito digitando o comando passwd na linha de comando, digitando a senha provisória, seguida da tecla Enter e da nova senha, duas vezes. Note que nada será exibido no terminal durante a digitação da nova senha!

    Recomendações:
    No Windows, utilize o Windows Terminal para linha de comando (disponível na Microsoft Store).
    No MacOS, utilize o Terminal ou iTerm.

    Caso precise de qualquer assistência adicional, não hesite em entrar em contato.

    Atenciosamente,
    Raphael Luiz Lobo da Silva Souza, Dr., TT-IV-A
    Laboratório de Estrutura e Evolução de Proteínas (LEEP/PSEL)
    Departamento de Microbiologia
    Instituto de Ciências Biomédicas
    Universidade de São Paulo
    Av. Prof. Lineu Prestes, 1374 - Ed. Biomédicas II - Sala 250 - 2º andar
    Tel: 3091-0891
    Cidade Universitária - CEP 05508-900 - São Paulo - SP - Brasil

    Se precisar de mais alguma coisa, estou à disposição!
    """

    try:
        server = smtplib.SMTP(smtp_server, smtp_port)
        server.starttls()  # Usar TLS (criptografia)
        server.login(from_email, from_password)

        msg = MIMEMultipart()
        msg['From'] = from_email
        msg['To'] = email
        msg['Subject'] = 'Acesso ao Servidor do CEPID B3'
        msg.attach(MIMEText(corpo_email, 'plain'))

        server.sendmail(from_email, email, msg.as_string())
        server.quit()
        print(f"E-mail enviado com sucesso para {email}")

        # Atualizar o status de "enviado" para "sim"
        update_sent_status(sheets_service, usuario)

    except Exception as e:
        print(f"Erro ao enviar o e-mail para {email}: {str(e)}")

# Função para atualizar o status de "enviado" para "sim" na planilha "Senhas"
def update_sent_status(sheets_service, usuario):
    spreadsheet_id = '1AYgGOApwYb7QsdNvAfXohGwFvIXaJzeFLol-eBFMFIM'  # ID do seu arquivo do Google Sheets
    range_to_check = "Senhas!A2:C"  # Intervalo na aba "Senhas" onde o login está na coluna A e o status na coluna C

    # Buscar todos os valores na coluna de logins e de status
    result = sheets_service.spreadsheets().values().get(spreadsheetId=spreadsheet_id, range=range_to_check).execute()
    values = result.get('values', [])

    if not values:
        print("Nenhum dado encontrado.")
        return

    # Procurar o login do usuário na coluna A (primeira coluna)
    for i, row in enumerate(values):
        if row[0] == usuario:
            # Encontramos o login, agora atualizamos a coluna de status para "sim"
            range_update = f"Senhas!C{i + 2}"  # A coluna de status está na coluna C, e as linhas começam de 2
            body = {
                "values": [["sim"]]
            }
            sheets_service.spreadsheets().values().update(
                spreadsheetId=spreadsheet_id,
                range=range_update,
                valueInputOption="RAW",
                body=body
            ).execute()
            print(f"Status de 'enviado' atualizado para 'sim' para o usuário {usuario}")
            break

# Função principal
def main():
    # Autenticar com o Google Drive e Sheets
    service, sheets_service = authenticate_google_drive()
    
    # O ID do arquivo Excel no Google Drive
    file_id = '1AYgGOApwYb7QsdNvAfXohGwFvIXaJzeFLol-eBFMFIM'  # ID do seu arquivo do Google Sheets
    
    # Baixar o arquivo Excel
    download_file(service, file_id)
    
    # Ler a segunda aba do arquivo Excel (.xlsx) baixado
    df = pd.read_excel('usuarios.xlsx', sheet_name="Enviar")  # Especificando a aba "Enviar"
    
    # Iterar sobre os dados
    for index, row in df.iterrows():
        nome = row['Nome']
        email = row['E-mail']
        senha = row['Senhas']
        usuario = row['Login']
        
        # Enviar o e-mail para cada usuário
        send_email(nome, email, senha, usuario, sheets_service)

    # Deletar os arquivos após a execução
    if os.path.exists("token.json"):
        os.remove("token.json")
        print("Arquivo token.json removido.")
    if os.path.exists("usuarios.xlsx"):
        os.remove("usuarios.xlsx")
        print("Arquivo usuarios.xlsx removido.")

if __name__ == "__main__":
    main()
